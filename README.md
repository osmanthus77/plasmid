# 质粒聚类
对质粒进行序列聚类，获得可聚质粒的共有序列，以便临床耐药性鉴定。   

## 1 NCBI获取RefSeq
```bash
mkdir ~/project/plasmid
cd ~/project/plasmid

rsync -avP ftp.ncbi.nlm.nih.gov::refseq/release/plasmid/ RefSeq/

gzip -dcf Refseq/*.genomic.gbff.gz > genomic.gbff
perl ~/project/plasmid/scripts/gb_taxon_locus.pl genomic.gbff > refseq_id_seq.csv
rm genomic.gbff

gzip -dcf RefSeq/plasmid.1.1.genomic.fna.gz | grep "^>" | head -n 5

faops n50 -S -C RefSeq/*.genomic.fan.gz

gzip -dcf RefSeq/*.genomic.fna.gz > RefSeq/plasmid.fa
```

## MinHash获取非冗余质粒non-redundant plasmids
```bash
mkdir ~/project/plasmid/nonredundant
cd ~/project/plasmid/nonredundant

faops size ../RefSeq/plasmid.fa > ../refseq.sizes

# 统计序列长度小于2000的序列数量
tsv-filter ../refseq.sizes --le 2:2000 | wc -l

```
