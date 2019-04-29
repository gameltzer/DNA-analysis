# Steps for assessment
* Gabriel Meltzer  
* Started 4/22/2019

---


1. ~~find ways to automate the downloading of the various links required.
	(urllib)~~
2. ~~Implement downloading.~~
3. ~~Illustrate method in Readme.~~
4. ~~understand the format of the data and document understanding~~
5. ~~~interpret instructions.~~~
6. ~~~Decide whether or not it should be placed into an SQLLITE database right away, or later. Decision: do it at the end:~~~
7. ~~~Assumption: work out how to extract information. Decision: pandas.~~~ 
8. ~~~Assumption: exclude entries listed in evidence file that do not relate to therapies.~~~
9. ~~~Assumption: exclude entries in variants files that are not missense variants. ~~~
10.~~~ Assumption: exclude entries in variants file that is combined missense and another type.~~~ 
11. ~~~Assumption: Exclude variatns that don't relate to a single AA position (some questions exist about this)~~~d
12.~~~ Assumption: exclude all evidence that does not relate to missense variants (might require joins)~~~
13. ~~~Assumption: remove all variants removed by the remaining evidence (might require joins)~~~
14. ~~~Assumption: create evidence_drug table ; one drug should exist per row; map to drugs column.~~~
15. ~~~Assumption: create variant_aliases in the tabl; map to variant_id. ~~~
16. ~~~extract amino acid position from each variant and make it a new field.~~~
17. ~~~ Output the DB into a file.~~~
18. Incorporate README contents into comments, especially concerning regular expressions. Revise comments and style of code, and fill in TODO statements relevant to comments.
19. Check database output to make sure that nothing is ommitted. 
20. revise overview.
21. Check regular expressions and other aspects of code, and add things specified in TODO statements. )
22. Rerun code to make sure none of the stylistic revisions or other additions created any errors.

(time as of Saturday is 30 hours)
