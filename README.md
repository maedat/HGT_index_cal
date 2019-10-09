
# HGT_index_cal; calculater for h and hA index.
====

HGT_index_cal calculte Horizontal gene transfer index "h" and the modified index "hA" from the blast result files.  The "h" index was calculated as the difference in bitscores between best prokaryote and best eukaryote matches in blast alignments, and the "hA" was obtained as the bitscores difference between best host- and best donar-database matches. 


## Requirement

 R (>3.6.1). (Tidyverse packages)

## Usage
```sh
Usage: aliva.py Rscript HGT_index.R dataA1 dataA2 dataB1 dataB2 type_name

dataA1-2; blast output (fmt6 format) for hA index calculation 
(e.g. A1=result against algal database, A2=result against animal database)

dataB1-2; blast output (fmt6 format) for h index calculation 
(e.g. B1=result against Prokaryote database, B2=result against Eukaryote database)

```

## Demo
```sh

Rscript HGT_index_cal.R  \
test/Apca_vs_dbAlgae.200.out  \
test/Apca_vs_dbAnimal.200.out   \
test/Apca_vs_dbProk.200.txt \
test/Apca_vs_dbEuk.200.txt \
Apca

```


## Licence
[MIT License](http://opensource.org/licenses/mit-license.php)

## Author
[Taro Maeda](https://github.com/maedat)
