# Skills Assessment (consider changing title and filename)
* Gabriel Meltzer
* 4/23/2019
* Python 2.7.15 (64 bit)
* Ubuntu 18.04.1 LTS and Windows 10 Home Version 1809 (see General Notes)
---

## General Notes
The software is being developed on a Lenovo Ideapad Flex 4 with Windows 10 Home Version 1809. To guarantee that it works across platforms, the code was also run on an Ubuntu 18.04.1 LTS Linux system through the use of an Oracle VirtualBox 5.2.8 r121009 ( Qt5.6.2) virtual machine.

## Downloading The Files
Downloading the files required for the assignment did not require anything beyond the standard library. The `urllib` package provided a function adequate for the task required. The function  `urlretrieve` allowed for a url to retrieve the data to be specified, as well as an argument specifying where it should be downloaded on the local file system (which includes the filename the file will be saved as). 

Since it was useful to also print status messages describing the progress of the operations, a new function `downloadTSV` was written. The argument `location` corresponded to the url sought, and the argument `name` corresponded to the filename (and location on the local system). This information was simply passed into the `urlretrieve` function. The value returned by the `urlretrieve` function was saved so that it could be used in the status message. 

## Extracting
The simplest, most straightforward thing to do would probably be to only load the data into SQLite once all the operations are performed. This will probably help reduce errors and better guarantee that only the desired information is the database. It is possible to just use text-processing to analyze the files, and then write all the transformations in intermediate files. Howver, the high number of operations that need to be performed means that keeping track of all the files will be tricky, and of course, if someone accidentally deletes a file while the operation is running, this will prevent the program from operating correctly. Thus, it would be better to do this within the memory of the program rather than saving to disk, which is also likely to prove faster. Thus, a decision was made to use the `pandas` library to perform these operations. 

The libary `pandas`  includes a function `read_csv` that can also be used to read tsv files. By changing the separator to a tab through the use of `\t`, data from tsv files can be easily extracted. The function returns a DataFrame that can then be more easily manipulated within python. For the purposes of this software, a function called `extract` was written which includes two calls to `read_csv` (one for each file). Some issues with the data required additional columns to be specified through the usecols argument for one of the files, so that the parser doesn't complain when more columns than expected are found. The two DataFrames representing the contents of the two files are then returned in a tuple. 

## Transforming

This is where the data is "shaped" to get it to the state it should be in for SQLite. The downloaded files are not themselves changed; rather, representations of the data are created in memory, and a series of operations are performed on them. Different sections will exist in this document discussing what was done for each of the steps. Note that all these operations are performed inside a `transform` function, which takes a tuple of two (with the data representations for evidence and variant), and also returns a tuple of two (with the transofrmed data represenations for evidence and variant)

### Transformations in Step 2 of instructions

The two operations that must be performed are to exclude all evidence that isn't related to drugs, and to exclude all variants that are not missense variations. There is a field in each data source that appears to correspond exactly to the information we are looking to use as a basis for filtering. 

For `evidence`, there is a field labled `drugs` which has the names of various drugs listed. The assumption is made that if there is evidence related to a drug for a particular entry, it will be listed here, with the name of the relevant drug(s) being listed in this column. It seem reasonable to assume that this column is the place where this information would be recorded. It is also assumed that entries that do not have a value in this column do not have any evidence pertaining to a particular drug, or it would be mentioned here. Thus, to exclude evidence that isn't related to drugs, all the rows are excluded where there is nothing available to indicate that any drug might be relevant for this particular row. Those are the rows, in the `drugs` column, with missing data, or `NaN`. `NaN` is treated as neither true nor false, so to exclude those rows, the `pandas` function `notna` is used to perform boolean indexing. The number of rows is now reduced from the original 2764 to 1817.

For `variants`, there are also missing `NaN` values in the DataFrame. These were filtered out since there is no way of knowing what they are; the assumption is that it is better to conclude that they are not missense variants, and thus be more conservative as to what is considered as a missense variant by requiring it to be stated explicitly.  The filtering out occurs with similar use of `notna` and Boolean indexing as described above. Since, for this step, the rows returned do not need to *only* be missense variants, the `str.contains` function of pandas was used rather than the equality operator. Thus, any rows that contain the phrase "missense_variant" are included in the DataFrame returned, and everything else is excluded. The assumptions was made that no misspellings or alternative descriptions of "missense_variant' occured. This operation reduced the number of rows from the original 754 as opposed to 2185. 

The variables `filteredEvidence` and `filteredVariant` that store the results of the operations are intended to be reused. Since the original data extracted is stored in separate `evidence` and `variant` DataFrames, a copy of the original is retained, and changes can be safely made since the orginal copy is still present to refer to. An alternative approach would be to create new variables for each operation, but it seems simpler to just resuse these variables once we perform the transformation rather than having to keep track of all the intermediate results. It also seems like a good to way insure that all the operations are cumulative, which is what is required. 

### Transformation in Step 3 of instructions

**Description of operations that must be performed**

For filtering the `filteredVariant` DataFrame so that we had missense variants and only missense variants, the operation is very similar to the one performed in the previous section, except, instead of `str.contains`, `str.match` is used, which is stricter and only includes items that are an exact match. Again, the assumption is made that there are no typos or alternate terminology for missense variants. The number of rows in the DataFrame, saved as `filteredVariant`, is 717.

Three things are worth investigating for the second part; hgvs expressions, allele reference id, and start position.

I encountered somewhat of a roadblock with the second operation, although I think a few assumptions might help me to to proceed. In any case, I have found that information about amino acid positons is found in two different columns, the `variant` column and the `hgvs_expressions` column. While I understand HGVS notation and how to interpret it, as well as information represented in a form like V600E, I noticed a lot of outliers that make things less straightforward.

The `variant` column in particularly is highly variable. Unable to find any specific clarification as to how non-specific variants would be represented at https://civicdb.org/help/variants/variants-overview, or through the HGVS website, my first approach was to just see if the variant column for a row contains a comma (since commas are often used in these fields to seperate multiple elements), using `str.contains`. It makes sense to assume that a non-specfic entry would have two or more elements sharing the same form as V600E (without repeating the same coordinates and amino acids, obviously), probably separated by a delimiter, and expecting that delimiter to be a comma seemed like a reasonable assumption.

Unfortunately, when I investigated a few of the results of the operation in the python interpreter, I did not find the results I had hoped to find, but I did see that a lot of the entries did not follow the pattern I had assumed. While sometimes the gene and additional information preceded the amino acid position information, this was not always the case.

Further investigation into the column revealed that information concerning amino acid positions was sometimes not even present. For instance, sometimes, the column would say something like "exon 2 mutation". Other times, a coordinate would be provided, but only the initial amino acid would be present, not the variant. 

#### Why is the information documented this way? 
One question I have would be concerning why this information was documented this way; as it seems like being more specific would be helfpul. For instance, if something is a missense variation, there should be an amino acid change, so it is unclear why this wouldn't be documented. 

While it might be possible to find or create a library to compute the amino acid change and display it, this would require information about the bases to be present, and that information is not always available in the reference_bases and variant_bases columns. The `hgvs_expressions` column also has the potential to help in that regard, but unfortunately, a lot of values are also missing in that column, so this is also not reliable. 

Would it be safe to just throw out this irregular information, or is there a benefit to keeping it? 

#### How are non-specific entries described, and what exactly do they mean? 

This was touched on before, but I was unable to find anything on either the CiViC website or http://www.hgvs.org/mutnomen/ expressing any standard for this kind of information. The chief question I would have to help resolve this issue would be to ask: How could I expect this information to be represented? 

I'm also curious as to the specific meaning of these non-specific entries... if these are multiple positions represented in the same row, why are they represented in the same row to begin with? Is it because, perhaps, they have the same effects?

#### A potential solution 

I think I may have found a way to proceed with this particular problem. I have opted to leave my questions in because I still happen to feel that this is not an optimal solution and I feel as though having additonal clarity with regards to the above question could perhaps assuage any doubts I might have. 

I opted to try using regular expressions with the `variant` field to see if  I could turn up anything useful. The function `str.contaions` allows the use of regular expressions if the boolean `regex` argument is set to True. The regex I used was "[ACDEFGHIKLMNPQRSTVWY]{1}[1234567890]+[ACDEFGHIKLMNPQRSTVWY]". The letters included in the character sets correspond to the one letter amino acid symbols,, and one and only one letter is expected for the amino acids at both ends of the "word". In between are of course the numbers indicating the actual position in the amino acid sequence. At least one number must be present here, or it would be impossible to know exactly which location this refers to. 

Accessing the variant field on the DataFrame that resulted from the reqular expression reveals that some of the rows have two codes containing amino acid positions separated by a **+**. Could this be the way the non-specific variants are represented? Given the absence of anything else that would appear to represent suhc a concept, I tentatively assume that this is the case, provided that the numerical portion is different for the different elements (which is not always the case).

Thus, I think a good approach to moving forward without having access to more information would be the following:

* DO NOT exclude rows that have a reference amino acid, followed by a sequence number, but without a variant amino acid. While I am unsure as to why this ifnormation would not be recorded, it does refer to a specific position without any ambiguity about the position. The only ambiguity is about the content of that position. Presumably, it would be more useful to include this information so that these positions can be investigated further before additional decisions are made. 

* Exclude rows that do not include data of the form "AAposAA" or "AApos" with the use of regular expressions, if the columns containing information related to the bases and the `hgvs_expressions` have missing data. In this case, it is impossible to relate it to a single AA position, so it seems reasonable to conclude that it is non-specific. At the very least, it is impossible to perform what is specified in step 7 with this data, so it seems reasonable to conclude that it is not relevant for our purposes.
* Exclude rows that contain a '+' sign in the variant column seperating two or more amino acid positions if the numerical portion of the position is not identical for the multiple positions. If this is the case, it is a nonspecific entry. An instance where the numerical portions are the same, however, don't strike me as being the same case; it seems somewhat analogous, in my mind, to having multiple alternate alleles at the nucleic acid level. It may, perhaps, even be a consequence of that. An additional assumption is made, then, that  something like  "V600E+V600M" should be regarded as a single AA position.

#####Implementation

The amino acid symbols were stored in a variable, `aminoAcidSymbols`. This variable was used, through string formatting, to create character classes for the regular expressions, as in `regexForAAPosition`. 

First, every entry that didn't contain amino acid positions in the `variant` field was filtered out through Boolean indexing and `str.contains`; `regexForAAPosition` was used as the regular expression, separated by a regex for the information in between the two amino acid positions. Raw strings are used because the regular expressions use backlashes. 

A regex for capturing the nonspecific amino acid positions, `nonspecificRegex` was created by combining two `regexForAAPosition` 

Then a regex is created to detect cases where there are two amino acid positions, but they are at the same coordinates ( like  "V600E+V600M"); recall that this is being treated as specific, rather than nonspecific. 

Finally, filtering occurs through the use of boolean indexing and `str.contains`. The codw within the brackets is actually in the form of an or statement, but care had to be taken to use the proper operators for pandas. Ultimately, the decision was made to include entries if the `variants` field either had two or more amino acid positions with the same corrdinates, or if it didn't contain two or more amino acids with different coordinates(i.e. it wasn't detected as nonspecific by the regex).

This should eliminate non-specific variant entries from the filtered results.

### Transformation in Step 4 of instructions

There are two operations involved. First, all evidence that does not relate to missense variants should be filtered out. Second, all variants not reference by the remaining evidence should be removed. 

The first operation requires 
*** enter info, and figure out how to verify. It is assumed that only the entries that remain in the `filteredVariant` DataFrame should be counted as missense variants, since the others were already removed from the data set and appear to not be of interest. Doing this will also guarantee that all the references in the `filteredEvidence` DataFrame are to entries in the `filteredVariant` DataFrame, avoiding potential confusion that might result in referring to something that was previously removed from the `filteredVariant` DataFrame. If it was necessary to include the items labeled as missense variants that were removed in previous operations, it would not be difficult to revise the code to accomodate this. The code should be identical, but it should be performed on the original `variant`DataFrame, rather than `filteredVariant`.


#### information to include that will be relevant later

** (relevant later)
.

**How can one make certain that the variant_id is included? ** why is it automatically removed? **

The other operation required in this particular step involves the amino acid position. There is no column directly corresponding to this. However, it appears that the `variant` column does contain this information, with amino acids represented by single letters, as in variant 5, where it is "C1156Y. The last "word" in the `variant` column contains two letters separated by a number. The first letter appears to describe the reference amino acid (what it would normally be), the number is the position in the sequence, and the last letter is the amino acid of the variant. 

 
Through looking at the `hgvs_expressions` for variant 5, the last element contains p-dot notation "NP_004295.2:p.Cys1156Tyr", which is more unambiguously an amino acid position; the *p* here stands for protein. It is found that it coresponds exactly to the suspected amino acid description mentioned before. The numbers are indentical, C corresponds to to Cysteine, and Y corresponds to Tyrosine. It is a reasonable to assume that this corespondence holds for all variants in this DataFrame, and that the word in the `variant` column is reliable place to obtain information about the amino acid position (in addition to containing this information in a more compact form than in the `hgvs_expressions` field)
**