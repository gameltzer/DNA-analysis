from urllib import urlretrieve
import pandas

def downloadTsv(location, name):
    print("downloading contents of {0}".format(location)) 
    result = urlretrieve(location,name)
    print("The file has been downloaded as {0}".format(result[0]))


evidenceUrl="https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv"
variantUrl="https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv"

evidenceFile="evidence.tsv"
variantFile="variant.tsv"

downloadTsv(evidenceUrl, evidenceFile)
downloadTsv(variantUrl,variantFile)

# This performs the extraction. A tuble is returned with the variant first, and the evidence second. 
def extract():
    # Some of the rows are read as having more columns than they expect. This issue appears to be resolved by including use cols, which apparently tells it to explicitly expect that greater number of columns. 
    # Only the python engine appears to know how to do this, which is why Python is specified. 
    evidence = pandas.read_csv(evidenceFile, "\t", engine="c")
    variant = pandas.read_csv(variantFile, "\t", engine="python", usecols=range(1,31))

    return (evidence,variant)    

evidence, variant = extract()
print(evidence.head())
print("n")
print(variant.head())
# print(variant.head() + "\n" )
# print(evidence.head())


# These are all the transform operations, including data cleaning. The evidence and the variant are returned as a tuple. Occasionally, missing values need to be handled, and operations are included to handle the
# missing values. The input ecpected is a tuple of two, with the first element containing the data of `evidence`,  and the second containing the data of `variant`. This is the same as return value of `extract()`. 
def transform(extractOutput):
    evidence, variant = extractOutput
    # This gets us only the evidence that concerns particular drugs. This is an assumption that is made. The variable filtered Evidence will be reused as we continue to transform and clean the data.
    filteredEvidence = evidence[evidence.drugs.notna()]
    # There are now 1817 rows in filteredEvidence as opposed to 2764
    # TODO: print out length and status of evidence at the moment. 
    
    # Since we have some NaN values present, we filter those out to avoid the possibility of screwing things up. 
    filteredVariant =  variant[variant.variant_types.notna()]
    # There are now 1309 rows in filtered variants as oppposed to 2185

    # This filters everything out that that is not a missense_varaint, although it may not only be a missense_variant. The assumption is made, of course, that there are no misspellings of missense_variant.
    filteredVariant = filteredVariant[filteredVariant.variant_types.str.contains("missense_variant")]
    # There are now 754 rows in filtered variants as opposed to 2185
    # TODO print length and status of variants at the moment. 

    # This filters the variants dataframe so that we only have rows where the variant_types consists of a missense and only a missense. 
    filteredVariant = filteredVariant[filteredVariant.variant_types.str.match("missense_variant")]
    # There are now 717 rows in filtered variants as opposeed t0 2185. 
    return (filteredEvidence, filteredVariant)
