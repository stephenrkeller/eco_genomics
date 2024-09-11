# Coding and data notes for the population genomics module

## Author: Steve Keller

### 09-10-2024: Intro to Centaurea GBS data and working with VCF files

We'll be analyzing the GBS data from 3 regions (EU, NE, PNW) starting today with variant call format files (VCF's)

-   We got to look at a fastq file using zcat \| head and discussed the Q-scores
-   Pick up with looking at fastp trim output on Th, and then proceed to mapping viewer (`tview`) using `load spack samtools`
-   Add student netIDs to the /netfiles/ecogen server before Thursday

### 09-11-2024: update to shared data access

* group management for //netfiles03/ecogen is handled by ETS, so we have to go through them to update enrollment via banner.  Decided this might not be most efficient way to proceed.
* Alternative is to copy dir's over to /gpfs1/cl/pbio3990 via VACC interactive terminal -- since both drives can be mounted and file/dir exchange can take place directly via `cp`
* So, moving forward, all paths will direct to /gpfs1/cl/pbio3990 for data access
