# practice

This program reads in a β-globin ([HBB, 11p15.4](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:4827)) cDNA sequence from a file and reports on the presence of the [rs334](https://www.ncbi.nlm.nih.gov/snp/rs334#clinical_significance) E7V variant. This variant produces a haemoglobin isoform known as HbS and causes [sickle cell anemia](https://www.orpha.net/consor/cgi-bin/Disease_Search.php?lng=EN&data_id=125&Disease_Disease_Search_diseaseGroup=232&Disease_Disease_Search_diseaseType=ORPHA&Disease(s)/group%20of%20diseases=Sickle-cell-anemia&title=Sickle%20cell%20anemia&search=Disease_Search_Simple) in affected individuals.
Coincidentally, this variant also confers resistance to malaria, which explains its relatively high prevalence in African populations ([MAF](https://en.wikipedia.org/wiki/Minor_allele_frequency): A=0.003-0.009).

In order to do that, the script transcribes the sequence into RNA, and translates the output in all 6 possible reading frames. 
Candidate amino-acid sequences featuring a stop codone are discarded.
The rest are compared against to targets: 
1) the wild-type HBB protein sequence, and 
2) the HBS isoform, harbouring `rs334`.

If there is a match with the HbS isoform, the script identifies the input as originating from a patient with sickle cell anemia.
