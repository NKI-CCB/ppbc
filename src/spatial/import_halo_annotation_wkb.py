import argparse
from dataclasses import dataclass, field
from datetime import datetime, timezone
from functools import reduce
from pathlib import Path
import re
from typing import List
import xml.etree.ElementTree as ET

import numpy as np
from shapely.geometry import LinearRing, Polygon
import shapely.ops
import shapely.validation
import shapely.wkb
import xarray as xr


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser("Convert HALO annotations to wkb")
    parser.add_argument("halo_xml", help="XML file exported from HALO")
    parser.add_argument("annotation_re",
            help="Regular expression to match the annotation to convert")
    parser.add_argument("out", help="WKB output file")
    args = parser.parse_args()
    args.halo_xml = Path(args.halo_xml)
    assert args.halo_xml.exists()
    assert args.halo_xml.is_file()
    args.out = Path(args.out)
    return args


def only_polygons(geometries):
    valid_geometries = list()
    for g in geometries:
        if g.area > 0:
            if isinstance(g, shapely.geometry.GeometryCollection):
                g = shapely.geometry.MultiPolygon([ig for ig in g.geoms
                                                   if isinstance(ig, shapely.geometry.Polygon)])
            valid_geometries.append(g)
    return valid_geometries


@dataclass
class Annotation:
    """A HALO Annotation, which is a named set of polygons."""

    name: str
    exterior_rings: List[LinearRing] = field(default_factory=list)
    interior_rings: List[LinearRing] = field(default_factory=list)
    holes: List[bool] = field(default_factory=list)

    def append(self, ring, hole):
        # Reverse ring direction if not concordant with hole status.  A hole should be
        # clockwise, while a ring with an external boundary should be anti-clockwise.
        if ring.is_ccw == hole:
            ring = LinearRing(ring.coords[::-1])
        if hole:
            self.interior_rings.append(ring)
        else:
            self.exterior_rings.append(ring)

    def as_multipolygon(self):
        # FIXME: This recreates only one level of holes, islands within holes are missed.
        # Merge polygons that are not holes. HALO output is not always a valid geometry,
        # so try to fix it.
        exterior_mpolygons = [shapely.validation.make_valid(Polygon(r))
                              for r in self.exterior_rings]
        interior_mpolygons = [shapely.validation.make_valid(Polygon(r.coords))
                              for r in self.interior_rings]
        # Remove lines and points, annotations should cover an area
        exterior_mpolygons = only_polygons(exterior_mpolygons)
        interior_mpolygons = only_polygons(interior_mpolygons)
        multi_polygon = reduce(lambda mp, p: mp.union(p), exterior_mpolygons)
        # Remove polygons that are holes
        multi_polygon = reduce(
            lambda mp, p: mp.difference(p), interior_mpolygons, multi_polygon
        )

        return multi_polygon


def str_to_bool(x):
    """Convert a string with 0 or 1 to a bool"""
    if x == "0":
        return False
    if x == "1":
        return True
    raise Excepion("Integer coded bool is not 0 or 1")


def read_annotation_xml(annot_xml):
    """Read a single annotation from a Annotation xml node"""
    annotation = Annotation(annot_xml.get("Name"))
    for region in annot_xml.find("Regions").iter("Region"):
        if region.get("HasEndcaps") != "0":
            raise Exception("Found HasEndcaps != 0")
        hole = str_to_bool(region.get("NegativeROA"))
        vertices = np.array(
            [(int(v.get("X")), int(v.get("Y"))) for v in region.find("Vertices").iter("V")]
        )
        if region.get("Type") == "Ellipse":
            if len(vertices) != 2:
                raise Exception("Ellipse should have two vertices")
            size_x = vertices[1][0] - vertices[0][0]
            size_y = vertices[1][1] - vertices[0][1]
            msize = (size_x + size_y) / 2
            mid_x = vertices[0][0] + (size_x / 2)
            mid_y = vertices[0][1] + (size_y / 2)
            circle = shapely.geometry.Point(mid_x, mid_y).buffer(msize/2)
            ellipse_poly = shapely.affinity.scale(circle, size_x / msize, size_y / msize)
            ring = ellipse_poly.exterior
        elif region.get("Type") == "Polygon":
            if len(vertices) < 3:
                raise Exception("Polygon should have at least three vertices")
            ring = LinearRing(vertices)
        else:
            raise Exception("Only polygons are supported")
        annotation.append(ring, hole)

    return annotation


def read_annotation_file(xml, annotation_re):
    """Read a HALO annotation from a xml file."""
    tree = ET.parse(xml)
    root = tree.getroot()

    if root.tag != "Annotations":
        raise Exception(f"{xml} is not an annotation file")

    annotation_re = re.compile(annotation_re)
    annotation = None
    for annot_xml in root.iter("Annotation"):
        if annotation_re.match(annot_xml.get("Name")):
            if annotation is not None:
                raise Exception("Multiple instances of annotation")
            annotation = read_annotation_xml(annot_xml)
    if annotation is None:
        raise Exception(f"Annotation {annotation} not found")
    return annotation


if __name__ == "__main__":
    args = parse_args()
    annotation = read_annotation_file(args.halo_xml, args.annotation_re)
    with open(args.out, "wb") as o:
        p = annotation.as_multipolygon()
        o.write(shapely.wkb.dumps(p))
