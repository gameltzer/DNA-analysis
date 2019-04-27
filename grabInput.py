from urllib import urlretrieve
import pandas

def downloadTsv(location, name):
    print("downloading contents of {0}".format(location)) 
    result = urlretrieve(location,name)
    print("The file has been downloaded as {0}".format(result[0]))

# These are useful for regular expressions. This is the single letter form. The lowercase characters are included just in case. 
aminoAcidSymbols= "ACDEFGHIKLMNPQRSTVWY"
aminoAcidSymbols= aminoAcidSymbols+aminoAcidSymbols.lower()

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
    # this allows us to use the index field, and rename it to variant_id. 
    variant = variant.reset_index()
    # This changes the column representing the old index back to variant_id.
    variant =variant.rename(columns={"index": "variant_id"})
    return (evidence,variant)    

    evidence, variant = extract()
    print(evidence.head())
    print("n")
    print(variant.head())
    # print(variant.head() + "\n" )
    # print(evidence.head())


    # These are all the transform operations, including data cleaning. The evidence and the variant are returned as a tuple. Occasionally, missing values need to be handled, and operations are included to handle the
    # missing values. The input ecpected is a tuple of two, with the first element containing the data of `evidence`,  and the second containing the data of `variant`. This is the same as return value of `extract()`. A tuple is returned consisting of the filtered evidence, the filtered variants, the evidence_drug DataFrame, and the variant_alias.
def transform(extractOutput):
    evidence, variant = extractOutput
        # This gets us only the evidence that concerns particular drugs. This is an assumption that is made. The variable filtered Evidence will be reused as we continue to transform and clean the data.
    filteredEvidence = evidence[evidence.drugs.notna()]
    # There are now va rows in filteredEvidence as opposed to 2764
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

    # This regex will be used to select all entries that have amino acid positional information in the variant field. 
    regexForAAPosition="[{0}]{{1}}[0-9]+[{0}]?".format(aminoAcidSymbols)
    filteredVariant = filteredVariant[(filteredVariant.variant.str.contains(pat=regexForAAPosition, regex=True)) | (filteredVariant.hgvs_expressions.notna())]
    
    # these next few lines make sure that only the specific variant entries remain.
    nonspecificRegex = regexForAAPosition+ r"\D+[\d\D]*" + regexForAAPosition
    specificRegexWCapture = r"[{0}]{{1}}([0-9]+)[{0}]?".format(aminoAcidSymbols)
    specificRegexWGroup = r"\D+[\d\D]*[{0}]{{1}}\1[{0}]?".format(aminoAcidSymbols)
    specificRegex = specificRegexWCapture +  specificRegexWGroup
    filteredVariant = filteredVariant[(filteredVariant.variant.str.contains(pat=specificRegex, regex=True)) | ~(filteredVariant.variant.str.contains(pat=nonspecificRegex, regex=True))  ]
   
   
    # TODO check to see if greediness or non-greediness is preferable
    #There are now 711 variants as opposed to 2185
    # TODO print length and status of variants at the moment.

   # This line makes sure filteredEvidence only contains the evidence that has variant ids that are in the variant dataframe table. Since all the variants in the variant dataframe are missense mutations, and sense they are the only missense variants we are interested in, this makes sure that only evidence pertaining to missense variants remains in the filteredEvidence dataframe. 
    filteredEvidence = filteredEvidence[filteredEvidence.variant_id.isin(filteredVariant.variant_id)]
     # TODO = write description of what the above does. 
    # TODO print length and status of evidence at the moment. 

    # TODO check to see if the filtering of evidence was correct. Are all the items in filteredEvidence also in filteredVariants? And did we check for which pertained to missense variants in the right way? 
    
    # This makes sure that only variants that are referred to in filteredEvidence are present in the dataframe. Since filteredEvidence refers to the same variant multile times, the size of filteredVariants will be smaller than filteredEvidence, instead of equal.
    filteredVariant = filteredVariant[filteredVariant.variant_id.isin(filteredEvidence.variant_id)]
    # TODO write description of what the above does.
    # TODO print length and status of variants at the moment.
    # TODO check to see that the filtering of variants was correct. We already know that evidence refers to the same variant many times, so the number will definitely not be 524. The number that we are getting is 204
    
    # Creates a new table with the drug for each evidence id (remember to convert to int when makeing it into sq) It is a DataFrame made from the evidence_id and drugs column of the filteredEvidecne DataFrame.
    evidence_drug = filteredEvidence.loc[:,["evidence_id","drugs"]]
    # evidence_drug.drugs.where(~(evidence_drug.drugs.str.contains(",")))

    # TODO wwrite description of what the above does. Include the assumption that the drugs will be delimited by commas and nothing else. There was a roadblock I encountered several times (an example of which can be seen in the row where evidence_id 305 is found.) which commas were within parentheses. It was unclear exactly what this means, and if only commas outside of parentheses should be counted. It looks like it may be some kind of typo, but perhaps it isn't. To resolve this, I would need to know the following: Why is the data in this particular format for those entries? Is it a mistake? Which parts represent drugs, and which parts  represent alternate names for the same drug? Before making a decision on the kind of regular expression to use to resolve this, it would be useful to know the answers to these questions. 

    # These DataFrames are the rows in evidence_drug with multiple drugs, and single drugs respectvely,
    #assuming that , is the delimiter for multiple entries shared bcross this field in the data set. They will be reassembled into a single DataFrame later.
    multiple_drugs = evidence_drug[evidence_drug.drugs.str.contains(",")]
    single_drugs = evidence_drug[~(evidence_drug.drugs.str.contains(","))]

    # This is an empty data frame for the loop to use. Later, it wil lbe assigned to the variable multiple_drugs. 
    temp_multiple = pandas.DataFrame()
    #  The field drugs is split around a , . New series are created for each drug, and these Series are appended to the data frame.
    for row in multiple_drugs.itertuples():
        drugsList =  row[2].split(",")
        for drug in drugsList: 
            rowToInsert = pandas.Series(name=row[0], data= [row[1], drug])
            temp_multiple = temp_multiple.append(rowToInsert)

    # This gives them the DataFrame the correct column names, and also makes sure the index is unique.  
    temp_multiple.columns =multiple_drugs.columns 
    temp_multiple = temp_multiple.reset_index() 

    # This reassigns the multiple_drugs variable to the DataFrame referred to by temp_multiple. The old, incorrect index is the first column, so we want to not include that. 
    multiple_drugs = temp_multiple.loc[:,["evidence_id", "drugs"]]

    # The multiple and single drugs DataFrames are reassembled into a single evidence_drug DataFrame. 
    evidence_drug =  multiple_drugs.append(single_drugs)
    # This makes sure we have an index with unique values, and this reoves the old, incorrect index. 
    evidence_drug = evidence_drug.reset_index()
    evidence_drug = evidence_drug.loc[:, ["evidence_id", "drugs"]]
    evidence_drug.rename(columns={"drugs":"drug"}).columns
    # TODO write detailed escriptoin of what was done
    # TODO print length and ststus of evidence_drugs
    # TODO Verify that the operation is working correctly. 
    
    # This creates the variant_alias from the DataFrame, and renames the variant_aliases column appropriately.
    variant_alias = filteredVariant.loc[:,["variant_id","variant_aliases"]]
    variant_alias = variant_alias.rename(columns={"variant_aliases":"alias"})


    # TODO: write description of what this does. Include the assumption that multiple variants aliases are delimited by a comma, and nothing else. The assumption is also made that those with missing values in the alias column can be removed from the DataFrame. Also, mention the potential roadblock of wondering if variant aliases are the same thing as amino acid positions, and if they are different, how are they different? The notation definitely appears to be amino acid positions in many cases.
    variant_alias = variant_alias[variant_alias.alias.notna()]
    multiple_aliases = variant_alias[variant_alias.alias.str.contains(",")] 
    single_alias = variant_alias[~(variant_alias.alias.str.contains(","))]

    temp_alias = pandas.DataFrame()
    for row in multiple_aliases.itertuples():
        aliasesList =  row[2].split(",")
        for alias in aliasesList: 
            rowToInsert = pandas.Series(name=row[0], data= [row[1], alias])
            temp_alias = temp_alias.append(rowToInsert)
    temp_alias.columns =multiple_aliases.columns 
    temp_alias = temp_alias.reset_index() 
    multiple_aliases = temp_alias.loc[:,["variant_id", "alias"]]
    variant_alias = multiple_aliases.append(single_alias)
    #TODO : print length and status of dataframe variant_alias
    #TODO : Confirm that this works.
    
    return (filteredEvidence, filteredVariant, evidence_drug, variant_alias)


# TODO write expression combining extract, transform and load so the statements are actually executed.