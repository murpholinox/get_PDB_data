---
title: "Obtención de datos del PDB"
layout: post
published: true
tag: PDB
---

## Objetivo

El propósito de este proyecto es comparar el daño por radiación en la proteína _x_ a un pH _i_ contra la misma proteína _x_ a un pH _i+n_; donde _n_ corre de uno hasta `MAX_VAL` (no determinado). Con esto en mente, el primer paso del análisis consiste en obtener una lista de proteínas que cristalicen en un intervalo de pH 'amplio'. Esta información se obtiene del PDB (del inglés _Protein Data Bank_), como se describe a continuación.

## Configuración

```{r, warning=FALSE, message=FALSE}
# Prepara el directorio de trabajo.
setwd("/home/murphy/Repos/doctorado/00_get_data")
# Carga el tidyverse.
library(tidyverse)
# Evalua o no bloques de código.
knitr::opts_chunk$set(eval = FALSE)
```

## Extracción de datos

### SIFTS

El proyecto [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/) (del inglés _structure integration with function, taxonomy and sequence_), provee, entre otras funciones, un mapeo entre entradas del PDB y UniProt. Dicho mapeo se usa para iniciar el análisis que se detalla a continuación.

```{bash, warning=FALSE, message=FALSE}
# Se descarga el archivo con el mapeo.
wget -O mappings.tsv.gz ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/uniprot_pdb.tsv.gz 
# Se descomprime el archivo.
gunzip mappings.tsv.gz
```

> El PDB maneja identificadores de cuatro caracteres. Por otra parte, UniProt usa identificadores de seis a diez caracteres (detalles [aquí](https://www.uniprot.org/help/accession_numbers)). El identificador en UniProt se denomina AC (del inglés _Accession Code_).

El archivo del mapeo `mappings.tsv` contiene dos columnas: el AC correspondiente a _x_ proteína y los identificadores del PDB pertenecientes a dicho AC, que a su vez corresponden a la proteína _x_. El siguiente bloque suma el número de entradas en el PDB para cada AC.

```{bash, warning=FALSE, message=FALSE}
# Elimina el cabezal.
tail -n +3 mappings.tsv > mappings_nohdr.tsv
# Cuenta entradas al contar ';'s.
awk '{print $0","gsub(/;/ ,"")}' mappings_nohdr.tsv> mappings_nohdr_counts.tsv
```

El archivo de salida `mappings_nohdr_counts.tsv` contiene ahora tres columnas; donde, a diferencia del anterior, se añadió la suma del número de entradas en el PDB de cada AC (por simplicididad, la `SUMA`). El siguiente bloque de código ayuda a ordenar los mapeos por `SUMA`.

```{r, warning=FALSE, message=FALSE}
uniprot_mappings <- read.csv("~/Repos/doctorado/00_get_data/mappings_nohdr_counts.tsv", header=FALSE)
# Ordena y obtiene los primeros cien.
top100_R <- head(arrange(uniprot_mappings, desc(uniprot_mappings$V2)), n=100)
# Los primeros seis caracteres de c1 es el AC. 
col_one <- str_sub(top100_R$V1, start = 1, end = 6)
# La `SUMA` es c2.
col_two <- top100_R$V2
top100 <- cbind(col_one, col_two)
write.csv(top100,"top100.csv", row.names = FALSE)
```

El archivo de salida `top100.csv` contiene dos columnas; el AC y la `SUMA` de las cien proteínas con mayor número de entradas en el PDB.

```{bash, warning=FALSE, message=FALSE}
# Elimina el cabezal.
tail -n +2 top100.csv | sed 's/"//g' > a
mv a top100_fx.csv
# Obtiene ACs.
awk -F "," '{print $1}' top100_fx.csv > top100_fx_uac.lst
# Convierte de lista a línea.
tr '\n' ' ' < top100_fx_uac.lst > top100_fx_uac.ln
```

El archivo de salida `top100_fx_uac.ln` del bloque anterior contiene una lista, en orden descendiente, de los ACs de las cien proteínas con mayor número de entradas en el PDB. Cabe señalar que hasta el momento los datos usados son solo identificadores del PDB y de UniProt. Por simplicidad es necesario conseguir los nombres recomendados por UniProt para estas cien proteínas.

### Nombres

El siguiente código invoca un programa en `perl`, el cual obtiene información de UniProt a partir de una lista de ACs.
> Este programa se encuentra [aquí](https://www.uniprot.org/help/api_batch_retrieval).

```{bash, warning=FALSE, message=FALSE}
# Descarga el programa.
wget -O get_info.pl https://raw.githubusercontent.com/murpholinox/usefulscripts/master/uniprot_batch_retrieval.pl
# Instala requisitos para correr el programa.
# sudo dnf install 'perl(LWP::UserAgent)' 
# sudo dnf install perl-LWP-Protocol-https
chmod u+x get_info.pl
# Corre el programa.
perl ./get_info.pl top100_fx_uac.ln > top100_wholeinfo.txt
# Obtiene identificadores.
egrep "^ID|^AC|^DE   RecN|^OS" top100_wholeinfo.txt > top100_ids_acs_fnames_org.txt
# Obtiene nombres.
grep "^DE" top100_ids_acs_fnames_org.txt | awk -F "=" '{print $2}' | sed 's/Full=//'g | sed 's/;$/,/g' > top100_fnames.lst
# Obtiene organismos.
grep "^OS" top100_ids_acs_fnames_org.txt | grep -v "HIV-1" | sed 's/^OS   //'g > top100_org.txt
# Algunos nombres contienen comas, esto provoca problemas.
# La siguiente línea halla las nombres problemáticos.
# awk -F , 'NF != 2 ' < top100_fnames.lst
# Y la siguiente línea arregla este problema.
sed -e 's/antigen, A al/antigen A al/g ; s/xide synthase, brain/xide synthase brain/g ; s/antigen, B al/antigen B al/g ; s/rylase, muscle form/rylase muscle form/g ; s/rases I, II, and III sub/rases I II and III sub/g ; s/xidase, mitoc/xidase mitoc/g' top100_fnames.lst > top100_fnames_fx.lst
# Pega los nombres de las proteínas con su AC, `SUMA` y organismo.
paste top100_fnames_fx.lst top100_fx.csv > a
paste -d, a top100_org.txt > top100_final.csv
rm a
# Se descartan algunas proteínas.
egrep -v "ribosomal protein|roteasome subunit|III subunit|polyprotein" top100_final.csv > top_final.csv
```

El archivo de salida `top_final.csv` contiene los nombres de las proteínas con su AC, `SUMA` y organismo al cual pertenecen. Cabe resaltar que del _top_ cien se eliminaron algunas entradas de proteínas cuya estructura cuaternaria es compleja (ribosoma, proteosoma) y otras donde a partir de la información genética se pueden codificar múltiples proteínas. El número de proteínas restantes disminuyó a 47.

### Consulta en el PDB

Con la lista de ACs restantes se hace una consulta en el PDB.

```{bash, warning=FALSE, message=FALSE}
# Obtiene ACs.
awk -F "," '{print $2}' top_final.csv > top_final_uacs.lst
# Convierte de lista a línea.
tr '\n' ' ' < top_final_uacs.lst > top_final_uacs.ln 
```

Se postea nuestra consulta en `xml` [aquí](http://www.rcsb.org/pdb/software/rest.do).

```
<orgPdbCompositeQuery version="1.0">
 <queryRefinement>
  <queryRefinementLevel>0</queryRefinementLevel>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>
    <description>Simple query for a list of Uniprot Accession IDs: </description>
    <accessionIdList>P00698 	P61769 	P00918 	P00720 	P11838 	P00760 	P24941 	P00734 	P56817 	P06746 	P04439 	O60885 	P42212 	P02766 	Q6PJP8 	O95696 	P02185 	P07900 	Q15596 	P03372 	P61823 	P01887 	P00644 	P69905 	P0AEX9 	P29476 	Q6B0I6 	P68871 	P01308 	P18031 	Q9UIF8 	P68431 	Q16539 	P19491 	P37231 	P00489 	Q15788 	P0CG48 	P01889 	P22629 	P61626 	P00533 	P00800 	P62805 	P04637 	P00431 	P0ABE7</accessionIdList>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>1</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.NumberOfEntitiesQuery</queryType>
    <description>Number of Entities Search : Entity Type=Protein Min Number of Entities=1 Max Number of Entities=1</description>
    <entity.type.>p</entity.type.>
    <struct_asym.numEntities.min>1</struct_asym.numEntities.min>
    <struct_asym.numEntities.max>1</struct_asym.numEntities.max>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>2</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
    <description>Experimental Method is X-RAY</description>
    <mvStructure.expMethod.value>X-RAY</mvStructure.expMethod.value>
    <mvStructure.expMethod.exclusive>y</mvStructure.expMethod.exclusive>
  </orgPdbQuery>
 </queryRefinement>
</orgPdbCompositeQuery>
```

Posteriormente, en la página de resultados, se seleccionan las siguientes columnas.

- Structure Summary
  - Resolution
  - Residue count
- Sequence 
  - Sequence
  - Chain length
  - External DB 
  - External DB Ref. ID
- Biological details
  - Source
  - Expression Host
  - EC number
- Sequence Clusters 
  - UniProt Accession Code 
  - UniProt Rec. Name
- Crystallization
  - pH value
  - Crystallization method
  - Crystal Growth Details
- Unit cell dimensions
  - Space group
-Primary citation
  - DOI
  - Title
    
Y el archivo resultante se almacena como `datos.csv`.








