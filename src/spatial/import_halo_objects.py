import csv
from dataclasses import dataclass
from datetime import datetime
from functools import partial
from pathlib import Path, PureWindowsPath
import re
import warnings

import click
import numpy as np
import netCDF4

cell_parts = ['nucleus', 'cytoplasm', 'membrane']

area_re = re.compile('((\w| )+) Area \((.*)\)')
dye_positive_cells_re = re.compile('(([a-zA-Z0-9]+ \(Opal [0-9]+\))|DAPI) Positive Cells')

class AutoIndexDict(dict):
    """A mapping of values to generated indices"""

    def __init__(self, *args, **kwargs):
        self.next_index = 0
        dict.__init__(self, *args, **kwargs)

    def __getitem__(self, key):
        val = self.get(key, None)
        if val is not None:
            return val
        else:
            val = self.next_index
            self.next_index += 1
            self[key] = val
            return val

def create_dimension_variable(g, name, values):
    """Make a netCDF dimension variable"""
    g.createDimension(name, len(values))
    g.createVariable(name, str, (name, ))
    for i, v in enumerate(values):
        g[name][i] = v

def create_string_variable(ds, name, values, dim):
    """Make a netCDF string variable"""
    ds.createVariable(name, str, (dim, ))
    for i, v in enumerate(values):
        ds[name][i] = v

@dataclass(frozen=True)
class Sample:
    t_number: str
    stain: str

    @classmethod
    def parse_filename(cls, filename):
        split = filename.split('_', 3)
        return cls(split[0], split[1])

def read_sample_summary(summary_fn):
    with open(summary_fn, encoding='utf-8-sig', newline='') as f:
        summary_reader = csv.DictReader(f)
        summary_header = summary_reader.fieldnames
        sample_summary = next(summary_reader)
        sample = Sample.parse_filename(sample_summary['Image Tag'])
        if next(summary_reader, None) is not None:
            raise Exception("More than one slide in a summary file")

    area_matches = (area_re.match(h) for h in summary_header)
    classes = sorted({m[1] for m in area_matches if m is not None}
                     - {'Classified', 'Avg Cell', 'Avg Nucleus', 'Avg Cytoplasm'})
    dye_matches = (dye_positive_cells_re.match(h) for h in summary_header)
    dyes = sorted({m[1] for m in dye_matches if m is not None})

    delimiter = ';'
    analysis_inputs = sample_summary['Analysis Inputs'].split(delimiter)
    analysis_inputs = (v.split(':', 1) for v in analysis_inputs)
    analysis_inputs = {v[0]: v[1] for v in analysis_inputs}
    analysis_inputs['Class List'] = analysis_inputs['Class List'].split(' | ')
    sample_summary['Analysis Inputs'] = analysis_inputs
    # Check that metadata in Analysis Inputs column matches column names
    negative_classes = {'Glass', 'Necrosis'}
    if set(analysis_inputs['Class List']) - negative_classes != set(classes) - negative_classes:
        raise Exception(
            "Classes do not match: [" +
            ", ".join(set(analysis_inputs['Class List'])) +
            "] versus ["
            ", ".join(set(classes)) +
            "]")

    return sample, sample_summary, dyes, classes

class Variable:
    def __init__(self, halo_name, nc_type, name=None, dims=('cell', ), units=None):
        self.halo_name = halo_name
        self.nc_type = nc_type
        if name is not None:
            self.name = name
        else:
            self.name = halo_name.lower().replace(' ', '_')
        self.dims = dims
        self.units = units
        self.fill_value = None

    def finalize(self, ds):
        if self.units is not None:
            ds[self.name].units = self.units

    def __str__(self):
        return f"Variable({self.halo_name}, nc_type={self.nc_type})"

def parse_uint(x, max_value = 2**32):
    x = int(x)
    assert x >= 0
    assert x < max_value
    return x

def parse_int(x, max_value = 2**31):
    x = int(x)
    assert x < max_value
    assert x > -max_value
    return x

class NumericVariable(Variable):
    
    def __init__(self, halo_name, nc_type, name=None, **kwargs):
        super().__init__(halo_name, nc_type, name, **kwargs)
        if nc_type in ['u1', 'u2', 'u4', 'u8']:
            self.fill_value = 2**(int(nc_type[1])*8)-1
            self.parser = partial(parse_uint, max_value=self.fill_value-1)
        elif nc_type in ['i1', 'i2', 'i4']:
            self.fill_value = 2**(int(nc_type[1])*8 - 1)
            self.parser = partial(parse_int, max_value=self.fill_value-1)
        elif nc_type in ['f4', 'f8']:
            self.parser = float
            self.fill_value = np.nan
        else:
            raise Exception(f"Unknown variable type {nc_type}")
    

class StringVariable(Variable):

    def __init__(self, halo_name, name=None, **kwargs):
        super().__init__(halo_name, str, name, **kwargs)
        self.parser = str

class FlagVariable(Variable):

    def __init__(self, halo_name, name=None, **kwargs):
        super().__init__(halo_name, 'u1', name=name, **kwargs)
        self.flag_indices = AutoIndexDict()
        self.parser = self.flag_indices.__getitem__

    def finalize(self, ds):
        super().finalize(ds)
        assert len(self.flag_indices) < 2**8
        ds[self.name].setncattr('flag_meanings', list(self.flag_indices.keys()))
        flag_values = np.array(list(self.flag_indices.values()), dtype='u1')
        ds[self.name].setncattr('flag_values', flag_values)

# Variables in the main objects file
variables = [
    NumericVariable('XMin', 'i4'),
    NumericVariable('XMax', 'i4'),
    NumericVariable('YMin', 'i4'),
    NumericVariable('YMax', 'i4'),
    NumericVariable('Object Id', 'u4'),
    NumericVariable('Positive Classification', 'u1', dims=('cell', 'dye')),
    NumericVariable('Cell Area (µm²)', 'f4', name='cell area', units='um2'),
    NumericVariable('Cytoplasm Area (µm²)', 'f4', name='cytoplasm area', units='um2'),
    NumericVariable('Nucleus Area (µm²)', 'f4', name='nucleus area', units='um2'),
    NumericVariable('Nucleus Perimeter (µm)', 'f4', name='nucleus perimeter', units='um'),
    NumericVariable('Nucleus Roundness', 'f4'),
    FlagVariable('Classifier Label'),
]
for cell_part in cell_parts:
    variables.append(NumericVariable(f"Positive {cell_part.title()} Classification", 'u1',
                     dims=('cell', 'dye')))
    variables.append(NumericVariable(f"{cell_part.title()} Intensity", 'f4',
                     dims=('cell', 'dye')))


def initialize_output(ds, sample, summary, variables, classes, dyes):
    n_cells = int(summary['Total Cells'])
    ds.createDimension('cell', n_cells)
    ds.createDimension('sample', 1)
    create_dimension_variable(ds, 'class', classes)
    ds.createDimension('dye', len(dyes))
    
    ds.createVariable('t_number', str, ('sample', ))
    ds['t_number'][0] = sample.t_number
    ds.createVariable('stain', str, ('sample', ))
    ds['stain'][0] = sample.stain
    ds.createVariable('image_location', str, ('sample', ))
    ds['image_location'][0] = summary['Image Location']
    ds.createVariable('image_tag', str, ('sample', ))
    ds['image_tag'][0] = summary['Image Tag']
    analysis_start_time = datetime.strptime(summary['Start Time'], '%d-%m-%Y %H:%M:%S')

    ds.analysis_start_time = analysis_start_time.isoformat()

    ds.createVariable('area', 'f8', ('class', 'sample'))

    ds['area'].units = 'mm2'
    for class_index, c in enumerate(classes):
        ds['area'][class_index] = summary[f'{c} Area (mm²)']

    create_string_variable(ds, 'dye', dyes, 'dye')
    
    for var in variables:
        if var.dims == ('cell',):
            chunksizes = (min(n_cells, 10000), )
        elif var.dims == ('cell', 'dye'):
            chunksizes = (min(n_cells, 10000), len(dyes))
        else:
            chunksizes = None
        ds.createVariable(var.name, var.nc_type, var.dims, chunksizes=chunksizes, zlib=True,
                          fill_value = var.fill_value)
    
    for k, v in summary['Analysis Inputs'].items():
        k = k.lower().replace(' ', '_')
        if k.startswith('dye'):
            continue  # FIXME: Dye variables are repeated
        ds.setncattr(k, v)



def write_rows(rows, start, ds, variables, dyes):
    for var in variables:
        if var.dims == ('cell',):
            ds[var.name][start:start+len(rows)] = [var.parser(r[var.halo_name]) for r in rows]
        elif var.dims == ('cell', 'dye'):
            for dye_idx, dye in enumerate(dyes):
                vn = f"{dye} {var.halo_name}"
                if  vn in rows[0]:
                    ds[var.name][start:start+len(rows), dye_idx] = [var.parser(r[vn]) for r in rows]
                # TODO: More efficiently encode missing dyes cell location combinations



@click.command()
@click.argument('summary_csv', type=click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.argument('halo_objects_csv', type=click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.argument('out', type=click.Path(exists=False, dir_okay=False, resolve_path=True))
def import_halo_objects(summary_csv, halo_objects_csv, out, sep=','):
    sample, summary, dyes, classes = read_sample_summary(summary_csv)

    buffer = []
    buffer_start = 0
    buffer_len = 1000
    with open(halo_objects_csv, encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)

        first_row = next(reader) 
        image_location = first_row['Image File Name']
        if PureWindowsPath(image_location).name != PureWindowsPath(summary['Image Location']).name:
            raise Exception("Summary file not about the same sample as the object file")

        ds = netCDF4.Dataset(out, 'w', format='NETCDF4')
        initialize_output(ds, sample, summary, variables, classes, dyes)
        buffer.append(first_row)
        for row in reader:
            assert row['Image File Name'] == image_location  # All data is about the same sample
            buffer.append(row)
            if len(buffer) % buffer_len == 0:
                write_rows(buffer, buffer_start, ds, variables, dyes)
                buffer_start += buffer_len
                buffer = []
        if len(buffer) > 0:
            write_rows(buffer, buffer_start, ds, variables, dyes)
        
        for var in variables:
            var.finalize(ds)

        ds.close()

if __name__ == '__main__':
    import_halo_objects()
