# Synteny_btw_regions

Series of commands that allow to analyse the synteny between two regions. It also allows to detect HGT if large Insertion/delitions between both regions.


# 1. Identify two regions to compare


- Extract region to compare from genome assembly of the first strain

      # create db 
      makeblastdb -dbtype nucl -in STRAIN_ASSEMBLY -parse_seqids -out STRAINdbFILE1_db 
      
      # extractregion with coordinates
      blastdbcmd -entry 'SCAFFOLD' -db STRAINdbFILE1_db -range minPosition-maxPosition > Destination/FileStrain1.fa


- Blast target region on second strains and extract sequence
      
      # create db 
      makeblastdb -dbtype nucl -in STRAIN_ASSEMBLY2 -parse_seqids -out STRAINdbFILE2_db
      
      # blast 1st region on the assembly of strain2
      blastn -db STRAINdbFILE2_db -outfmt 6 -evalue 1e-8 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out blastStrain1_in_Strain2.xml -query Destination/FileStrain1.fa
      
      # extractregion with coordinates
      blastdbcmd -entry 'SCAFFOLD' -db STRAINdbFILE2_db -range minPosition-maxPosition > Destination/File2.fa
      
      
# 2.  Annotate region of both strains
            
- With prokka for prokaryotic genomes
 
            ~/software/prokka-1.12/prokka/bin/prokka --outdir Annotation_FOLDER --genus GENUS --species SPECIES --strain STRAIN --locustag GS_STRAIN --prefix FILEPREFIX_Prokka --rfam --usegenus Destination/File.fa

- For eukaryotes
      
     - Use genemark to detect Open reading frame (ORF) http://opal.biology.gatech.edu/GeneMark/
      
     - Use Blast2GO for annotate the ORF. (It is time consumming!!) https://www.blast2go.com/
      
# 3. Sinteny visualization       
            
