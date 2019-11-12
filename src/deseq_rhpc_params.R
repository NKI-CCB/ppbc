
#Due to issues to BLAS on the RHPC server, DESeq will parallelize itself even when no parallel backend is registered
#This behavior drastically reduces the speed of calls to DESeq
#To prevent this behavior, an environment variable needs to be set BEFORE R is started
#This eliminates the use of server-based RStudio

#The workaround is as follows:
#Open the bash terminal and log into harris/darwin
#Run this command:
#export OMP_NUM_THREADS=1
#Then run the script from the console, using Rscript yourscript.R

#You might be tempted to try this either from within a script or within RStudio
#Sys.setenv(OMP_NUM_THREADS=1)
#But it won't work because the environment variable needs to be set before R starts

