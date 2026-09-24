#script copied from https://github.com/cguccione/human_host_filtration/
# Feb 6th, 2025
#!/bin/bash

# Build human minimap2 databases
echo "Building human minimap2 databases"
THREADS="${1:-12}"

minimap2 -ax sr -t "${THREADS}" -d ref/human-GRC-db.mmi ref/GRCh38_latest_genomic.fna
minimap2 -ax sr -t "${THREADS}" -d ref/human-GCA-phix-db.mmi ref/human-GCA-phix.fna

# Remove large unneeded files
rm ref/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna ref/GRCh38_latest_genomic.fna ref/human-GCA-phix.fna

echo "Building pangenome minimap databases"
directory_path="ref/pangenomes"

for file in "$directory_path"/*
do
    if [ -f "$file" ]; then
        echo "Indexing $file"
        filename=$(basename "$file")
        mmi_name="${filename%.*}"
        minimap2 -d "$directory_path/$mmi_name.mmi" "$file"
    fi
done

echo "Done indexing"
