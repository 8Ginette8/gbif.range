# Download worldwide the observations of Panthera tigris and Ailuropoda melanoleuca
obs.pt = get_gbif("Panthera tigris")
obs.am = get_gbif("Ailuropoda melanoleuca")
\dontrun{
# Retrieve DOI for only one get_gbif() output
get_doi(obs.pt,title="GBIF_test1",description="A small example 1",
   source_url="https://example.com/",user="",pwd="") # Use your own GBIF credentials here

# Retrieve DOIs for several get_gbif() outputs
get_doi(list(obs.pt,obs.am),title="GBIF_test2",description="A small example 2",
   source_url="https://example.com/",user="",pwd="") # Use your own GBIF credentials here
}