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

faops n50 -S -C RefSeq/*.genomic.fna.gz

gzip -dcf RefSeq/*.genomic.fna.gz > RefSeq/plasmid.fa
```

## 2 MinHash获取非冗余质粒non-redundant plasmids
```bash
mkdir ~/project/plasmid/nonredundant
cd ~/project/plasmid/nonredundant

faops size ../RefSeq/plasmid.fa > refseq.sizes

# 统计序列长度小于2000的序列数量
tsv-filter ../refseq.sizes --le 2:2000 | wc -l

# 提取长度大于2000的序列
faops some ../RefSeq/plasmid.fa <(tsv-filter refseq.sizes --gt 2:2000) refseq.fa

# 得到sketch/构建本地数据库
mash sketch -k 21 -s 1000 -i -p 8 refseq.fa -o refseq.plasmid.k21s1000.msh

# split
mkdir -p job
faops size refseq.fa | cut -f 1 | gsplit -l 1000 -a 3 -d - job/

# 每个小文件分别生成mash sketch
find job -maxdepth 1 -type f -name "[0-9]??" | sort | 
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        faops some refseq.fa {} stdout | 
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
    '

# 计算每个小文件与refseq的Mash距离
find job -maxdepth 1 -type f -name "[0-9]??" | sort | 
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        mash dist -p 6 {}.msh refseq.plasmid.k21s1000.msh > {}.tsv
    '

# 提取距离小于0.01的部分(相似度高)
find job -maxdepth 1 -type f -name "[0-9]??" | sort | 
    parallel -j 16 '
        cat {}.tsv | 
            tsv-filter --ff-str-ne 1:2 --le 3:0.01
    ' \ 
    > redundant.tsv

head -n 5 redundant.tsv
cat redundant.tsv | wc -l

# 划分连通组件，同一连通组件内相似度高
cat redundant.tsv | 
    perl -nla -F"\t" -MGraph::Undirected -e '
        BEGIN {
            our $g = Graph::Undirected->new; 
        }
        
        $g->add_edge($F[0], $F[1]); 
        
        END{
            for my $cc ( $g->connected_components ){
                print join qq{\t}, sort @{$cc};
            }
        }
    ' \
    > connected_components.tsv

cat connected_components.tsv | perl -nla -F"\t" -e 'printf qq{%s\n}, $_ for @F' > components.list

wc -l connected_components.tsv components.list

faops some -i refseq.fa components.list stdout > refseq.nr.fa
faops some refseq.fa <(cut -f 1 connected_components.tsv) stdout >> refseq.nr.fa

rm -rf job
```

## 3 MinHash聚类
```bash
mkdir ~/project/plasmid/grouping
cd ~/project/plasmid/grouping

cat ../nonredundant/refseq.nr.fa | mash sketch -k 21 -s 1000 -i -p 8 - -o refseq.nr.k21s1000.msh

# 分割为小文件
mkdir -p job
faops size ../nonredundant/refseq.nr.fa | cut -f 1 | gsplit -l 1000 -a 3 -d - job/

# 得到nr的sketch
find job -maxdepth 1 -type f -name "[0-9]??" | sort | 
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}
        faops some ../nonredundant/refseq.nr.fa {} stdout | 
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
    '

# 计算mash距离
find job -maxdepth 1 -type f -name "[0-9]??" | sort | 
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        mash dist -p 6 {}.msh refseq.nr.k21s1000.msh > {}.tsv
    '

# 小文件的距离合并为一个文件
find job -maxdepth 1 -type f -name "[0-9]??" | sort | 
    parallel -j 1 '
        cat {}.tsv
    ' \ 
    > dist_full.tsv

# distance < 0.05
cat dist_full.tsv | tsv-filter --ff-str-ne 1:2 --le 3:0.05 > connected.tsv

head -n 5 connected.tsv
cat connected.tsv | wc -l

# 划分连通组件并按照大小分类保存
mkdir -p group
cat connected.tsv | 
    perl -nla -F"\t" -MGraph:Undirected -MPath::Tine -e '
        BEGIN {
            our $g = Graph::Undirected->new;
        }
        
        $g->add_edge($F[0], $F[1]);
        
        END {
            my @rare;
            my $serial = 1;
            my @ccs = $g->connected_components;
            @ccs = map { $_->[0] }
                sort { $b->[1] <=> $a->[1] }
                map { [ $_, scalar (@{$_} ) ] } @ccs
            for my $cc ( @ccs ) {
                my $count = scalar @{$cc};
                if ($count < 50 ) {
                    push @rare, @{$cc};
                }
                else {
                    path(qq{group/$serial.lst})->spew(map {qq{$_\n}} @{$cc});
                    $serial++;
                }
            }
            path(qq{group/00.lst})->spew(map {qq{$_\n}} @rare);
            
            path(qq{grouped.lst})->spew(map {qq{$_\n}} $g->vertices);
        }
    '







```