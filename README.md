# Synteny_btw_regions

Series of commands that allow to analyse the synteny between two regions. It also allows to detect HGT if large Insertion/delitions between both regions.


# 1. Identify two regions to compare

- Blast target region on second strains.


- Extract region to compare from genome assembly

      # create db 
      makeblastdb -dbtype nucl -in STRAIN_ASSEMBLY -parse_seqids -out STRAINdbFILE_db 
      
      # extractregion with coordinates
      blastdbcmd -entry 'SCAFFOLD' -db STRAINdbFILE_db -range minPosition-maxPosition > Destination/File.fa
      
      
- Annotate region

