-- mysql -u cod_reader -h sql.crystallography.net cod < count-COD-structures.sql
SELECT file,formula,svnrevision FROM data
WHERE
-- For sure we do not want retracted and erroneous structures. The
-- COD numbers however for these structures stay in the COD , to
-- prevent their repeated deposition. We filter them out:
(STATUS IS NULL OR
(status NOT LIKE "%retracted%" AND status NOT LIKE "%error%")
) AND
-- In this SELECT , we search for oxides (element 'O' in the
-- formula ):
formula REGEXP " F[0-9.]* " AND 
-- NOT(formula REGEXP "\\.") AND
NOT(
-- The excluded elements from the manuscript:
formula regexp " (N|O|P|S|Cl|As|Se|Br|Sb|Te|I|Bi|Po|At|C|H)[0-9.]* "
) AND
-- Some small number of structures landed in the COD several times;
-- when spotted such entries are marked as 'duplicates ' and can be
-- easily filtered away:
duplicateof IS NULL AND
-- Experiment conditions as mentioned by the Authors. Not that in
-- many cases temperatures and pressures for normal conditions are
-- not given , so we take the NULL structures as well:
(celltemp IS NULL OR celltemp BETWEEN 270 AND 370) AND
(diffrtemp IS NULL OR diffrtemp BETWEEN 270 AND 370) AND
(cellpressure IS NULL OR cellpressure < 1000) AND
(diffrpressure IS NULL OR diffrpressure < 1000) AND
-- For the purpose of the investigation , we need coordinates:
flags LIKE "%has coordinates%" AND
-- Some software can not deal with disorder. If this is the case we
-- can exclude such structures as well. In any case , their number
-- is nor large here:
flags NOT LIKE "%has disorder%" AND
-- We want only experimental structures , so let 's remove all
-- entries marked as "theoretical ":
(method IS NULL OR method != "theoretical")
-- AND
-- To reproduce this query tomorrow , we specify a particular
-- revision , since the COD is updated daily:
-- svnrevision <= 271589
;
