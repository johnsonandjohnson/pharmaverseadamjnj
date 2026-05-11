b.	PARAM = “End of Trt Status for study agent 2”
c.	AVALC = randomly assign values of “ONGOING”, “DISCONTINUED”, “COMPLETED”
d.	DSSCAT = “Study Agent 1”
5.	Create a record for each subject containing a record with PARAMCD=EOTS2STT and AVALC=”DISCONTINUED” with the following (note: subjects that do not have a record with PARAMCD=EOTS2STT and AVALC=”DISCONTINUED”  will not have record created):
a.	PARAMCD = “DCTS2RS”
b.	PARAM = “Reason for Treatment Discontinuation for study agent 2”
c.	AVALC = “Other”
d.	AVALCSP = “specify text”
e.	ASTDT = adsl.DCTDT
f.	ASTDTC = adsl.DCTDT converted to iso8601 format (e.g. yyyy-mm-dd; 2026-05-07) as a character variable
g.	ASTDY = adsl.DCTADY
6.	Create a record for each subject containing non-missing adsl.LTVISIT with the following (note: subjects that are missing adsl.LTVISIT will not have record created):
a.	PARAMCD = “LTVISTS1”
b.	PARAM = “Last Treatment Visit for study agent 1”
c.	AVALC = adsl.LTVISIT
7.	Create a record for each subject containing non-missing adsl.LTVISIT with the following (note: subjects that are missing adsl.LTVISIT will not have record created):
a.	PARAMCD = “LTVISTS2”
b.	PARAM = “Last Treatment Visit for study agent 2”
c.	AVALC = adsl.LTVISIT
8.	Create subject level variables
a.	S1EDT = adsl.TRTEDT
b.	S1EDTM = adsl.TRTEDTM
c.	S1EDY = adsl.TRTEDY
d.	S2EDT = adsl.TRTEDT
e.	S2EDTM = adsl.TRTEDTM
f.	S2EDY = adsl.TRTEDY
g.	S1SDT = adsl.TRTSDT
h.	S1SDTM = adsl.TRTSDTM
i.	S1SDY = adsl.TRTSDY
j.	S2SDT = adsl.TRTSDT
k.	S2SDTM = adsl.TRTSDTM
l.	S2SDY = adsl.TRTSDY
9.	Retrieve the following subject level variables from pharmaverseadamjnj and join them with the dataset after step #8 so that they appear on each record
a.	PPROTFL, SAFFL, RANDFL, SCRNFL, SCRFFL, RESCRNFL, FASFL, ENRLFL, AGE, AGEU, SEX, RACE, COUNTRY, RFICDT, SITEID, SUBJID, TRT01P, TRT01PN, TRT01A, TRT01AN, 