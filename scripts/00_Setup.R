# Download the data
# Data Source: https://www.nature.com/articles/s41598-020-64929-x
download.file(
    url = "https://single-cell-transcriptomics.s3.eu-central-1.amazonaws.com/course_data.tar.gz",
    destfile = "data.tar.gz"
)

# untar
untar("data.tar.gz")

# remove compressed files
file.remove("data.tar.gz")
