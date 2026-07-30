## Test environments

* local: Windows 10 x64, R 4.3.1

## R CMD check results

0 errors | 0 warnings | 4 notes

### Note 1: CRAN incoming feasibility

This is the first submission of gbif.range to CRAN.

The same note reports five URLs as possibly invalid, all returning HTTP 403:

    https://www.gbif.org
    https://www.gbif.org/developer/occurrence
    https://www.gbif.org/developer/summary
    https://www.iucnredlist.org
    https://www.iucnredlist.org/resources/spatial-data-download

All five have been checked manually and are valid and reachable in a browser.
Both gbif.org and iucnredlist.org return 403 to automated requests that do not
present a browser user agent. GBIF is the data source this package is built
around, so these references are essential and have no alternative canonical
form.

### Note 2: Installed size

    installed size is 5.6Mb
    sub-directories of 1Mb or more:
      doc       1.7Mb
      extdata   2.2Mb
      help      1.2Mb

The extdata directory holds the small spatial and occurrence files used by the
vignettes and examples. These allow a substantial part of the documentation to
run without contacting a remote server, which we consider worth the size. The
doc directory holds five vignettes covering the package workflows.

### Note 3: Unable to verify current time

This appears only on our local machine, which cannot reach the time server used
by the check. We do not expect it on the CRAN build machines.

### Note 4: Examples with elapsed time > 5s

The package retrieves occurrence records from GBIF and ecoregion layers from
external providers. For most of the examples listed, elapsed time is dominated
by waiting on network responses rather than by computation: get_status uses
0.5s of CPU for 13.6s elapsed, and obs_filter 6.8s of CPU for 38.1s elapsed.
All network-dependent examples are wrapped in \donttest{}.

In line with the CRAN policy on internet resources, get_gbif() retries and then
returns an empty result with an informative message when GBIF is unreachable,
rather than raising an error, and read_ecoreg() returns NULL with a message
when an ecoregion layer cannot be downloaded. The affected examples test their
inputs before proceeding, so a temporary outage at an external provider cannot
turn into a check error.

One example, get_doi(), is wrapped in \dontrun{} because it requires personal
GBIF account credentials and cannot be executed without them.
