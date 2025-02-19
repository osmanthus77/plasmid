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
MinHash原理：[通俗易懂理解Minhash算法](https://blog.csdn.net/a61022706/article/details/141640451)
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

# 计算每个小文件与refseq的 Mash距离
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
    ' > redundant.tsv

head -n 5 redundant.tsv
cat redundant.tsv | wc -l

# 构建无冗余网络（划分连通组件，同一连通组件内相似度高）
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
    ' > connected_components.tsv

# 列转行
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
        echo >&2 "==> {}"
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
    ' > dist_full.tsv

# distance < 0.05
cat dist_full.tsv | tsv-filter --ff-str-ne 1:2 --le 3:0.05 > connected.tsv

head -n 5 connected.tsv
cat connected.tsv | wc -l

# 划分连通组件并按照大小分组保存序列ID
mkdir -p group
cat connected.tsv | 
    perl -nla -F"\t" -MGraph::Undirected -MPath::Tiny -e '
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
                map { [ $_, scalar( @{$_} ) ] } @ccs;
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

# 未成组的
faops some -i ../nonredundant/refseq.nr.fa grouped.lst stdout | 
    faops size stdin | 
    cut -f 1 > group/lonely.lst

wc -l group/*

# 每个组件内比较相似性
find group -maxdepth 1 -type f -name "[0-9]*.lst" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        
        faops some ../nonredundant/refseq.nr.fa {} stdout | 
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
            
        mash dist -p 6 {}.msh {}.msh > {}.tsv
    '

# 每个连通组件内部层次聚类
find group -maxdepth 1 -type f -name "[0-9]*.lst.tsv" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"

        cat {} |
            tsv-select -f 1-3 |
            Rscript -e '\''
                library(readr);
                library(tidyr);
                library(ape);
                pair_dist <- read_tsv(file("stdin"), col_names=F);
                tmp <- pair_dist %>%
                    pivot_wider( names_from = X2, values_from = X3, values_fill = list(X3 = 1.0) )
                tmp <- as.matrix(tmp)
                mat <- tmp[,-1]
                rownames(mat) <- tmp[,1]

                dist_mat <- as.dist(mat)
                clusters <- hclust(dist_mat, method = "ward.D2")
                tree <- as.phylo(clusters)
                write.tree(phy=tree, file="{.}.tree.nwk")

                group <- cutree(clusters, h=0.2) # k=3
                groups <- as.data.frame(group)
                groups$ids <- rownames(groups)
                rownames(groups) <- NULL
                groups <- groups[order(groups$group), ]
                write_tsv(groups, "{.}.groups.tsv")
            '\''
    '

# subgroup
mkdir -p subgroup
cp group/lonely.lst subgroup/

find group -name "*.groups.tsv" | sort | 
    parallel -j 1 -k '
        cat {} | sed -e "1d" | xargs -I\[\] echo "{/.}_[]"
    ' |
    sed -e 's/\.lst\.groups_/_/' |
    perl -na -MPath::Tiny -e '
        path(qq{subgroup/$F[0].lst})->append(qq{$F[1]\n});
    '
 
# ignore small subgroups
find subgroup -name "*.lst" | sort |
    parallel -j 1 -k '
        lines=$(cat {} | wc -l)
        
        if (( lines < 5 )); then
            echo -e "{}\t$lines"
            cat {} >> subgroup/lonely.lst
            rm {}
        fi
    '

# 将前一步骤非冗余中的connected_components中内容添加到subgroup中
cat ../nonredundant/connected_components.tsv | 
    parallel -j 1 --colsep "\t" '
        file=$(rg -F -l  "{1}" subgroup)
        echo {} | tr "[:blank:]" "\n" >> ${file}
    '

# 删除重复
find subgroup -name "*.lst" | sort |
    parallel -j 1 '
        cat {} | sort | uniq > tmp.lst
        mv tmp.lst {}
    '

wc -l subgroup/* | sort -nr | head -n 100

wc -l subgroup/* | perl -pe 's/^\s+//' | tsv-filter -d " " --le 1:10 | wc -l

# 筛选subgroup中多余50条序列的组
wc -l subgroup/* | perl -pe 's/^\s+//' | tsv-filter -d " " --ge 1:50 | tsv-filter -d " " --regex '2:\d+' | sort -nr > next.tsv
wc -l next.tsv

rm -rf job
```

## 4 Plasmid：prepare预处理
### 拆分序列
```bash
mkdir ~/project/plasmid/GENOMES
mkdir ~/project/plasmid/taxon
cd ~/project/plasmid/grouping

echo -e "#Serial\tGroup\tCount\tTarget" > ../taxon/group_target.tsv

# 统计分组信息（编号、subgroup组名、序列数量、代表序列ID）
# 并按照subgroup分组顺序从refseq提取序列
cat next.tsv | cut -d " " -f 2 | 
    parallel -j 4 -k --line-buffer '
        echo >&2 "==> {}"
        
        GROUP_NAME={/.}
        TARGET_NAME=$(head -n 1 {} | perl -pe "s/\.\d+//g")
        
        SERIAL={#}
        COUNT=$(cat {} | wc -l)
        
        echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}" >> ../taxon/group_target.tsv
        
        faops order ../nonredundant/refseq.fa {} stdout | 
            faops filter -s stdin stdout \
            > ../GENOMES/${GROUP_NAME}.fa
    '

# 统计subgroup中序列大小
cat next.tsv | cut -d" " -f 2 |
    parallel -j 4 -k --line-buffer '
        echo >&2 "==> {}"
        GROUP_NAME={/.}
        faops size ../GENOMES/${GROUP_NAME}.fa > ../taxon/${GROUP_NAME}.sizes
    '

# optional：RepeatMasker
#egaz repeatmasker -p 16 ../GENOMES/*.fa -o ../GENOMES/

# 按照序列ID拆分序列split-name
find ../GENOMES -maxdepth 1 -type f -name "*.fa" | sort |
    parallel -j 4 '
        faops split-name {} {.}
    '

# 拆分后序列，一条序列单独一个文件夹
find ../GENOMES -maxdepth 2 -mindepth 2 -type f -name "*.fa" | sort |
    parallel -j 4 '
        mkdir -p {.}
        mv {} {.}
    '
```

### prepseq预处理（preparing steps for lastz）
```bash
cd ~/project/plasmid

cat taxon/group_target.tsv | sed -e '1d' |
    parallel --colsep '\t' -no-run-if-empty --line-buffer -k -j 4 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"
        
        for name in $(cat taxon/{2}.sizes | cut -f 1); do
            egaz prepseq GENOMES/{2}/${name}
        done
    '
```

### 筛选序列长度异常值
```bash
cd ~/project/plasmid

# 统计next.tsv序列数量、长度之和
cat taxon/*.sizes | cut -f 1 | wc -l
# cat taxon/*.sizes | cut -f 2 | paste -sd+ | bc
awk '{sum+=$2} END {print sum}' taxon/*.sizes

cat taxon/group_target.tsv | sed -e '1d' |
    parallel --colsep '\t' --no-run-if-empty --line-buffer -k -j 4 '
        echo -e "==> Group: [{2}]\t Target: [{4}]"
        
        median=$(cat taxon/{2}.sizes | datamash median 2)
        mad=$(cat taxon/{2}.sizes | datamash mad 2)
        lower_limit=$( bc <<< " (${median} - 2 * ${mad}) / 2" )

#        echo $median $mad $lower_limit
        lines=$(tsv-filter taxon/{2}.sizes --le "2:${lower_limit}" | wc -l)
        
        if (( lines >0 )); then
            echo >&2 "     $lines lines to be filtered"
            tsv-join taxon/{2}.sizes -e -f <(
                    tsv-filter taxon/{2}.sizes --le "2:${lower_limit}"
                ) \
                > taxon/{2}.filtered.sizes
            mv taxon/{2}.filtered.sizes taxon/{2}.sizes
        fi
    '

cat taxon/*.sizes | cut -f 1 | wc -l
# cat taxon/*.sizes | cut -f 2 | paste -sd+ | bc
awk '{sum+=$2} END {print sum}' taxon/*.sizes
```

### Rsync to hpcc
```bash
rsync -avP ~/project/plasmid/ wangq@202.119.37.251:data/plasmid
```

## Plasmid:run
```bash
cd ~/project/plasmid


cat taxon/group_target.tsv | sed -e '1d' | grep "^53" |
    parallel --colsep '\t' --no-run-if-empty --line-buffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget:[{4}]\n"
        
        egaz template GENOMES/{2}/{4} \
            $(cat taxon/{2}.sizes | cut -f 1 | grep -v -x "{4}" | xargs -I \[\] echo "GENOMES/{2}/[]") \
            --multi -o groups/{2}/ \
            --order \
            --parallel 24 -v

#           bash groups/{2}/1_pair.sh
#           bash groups/{2}/3_multi.sh

            bsub -q -mpi -n 24 -J "{2}-1_pair" "bash groups/{2}/1_pair.sh"
            bsub -w "ended({2}-1_pair)" \
                -q mpi -n 24 -J "{2}-3_multi" "bash groups/{2}/3_multi.sh"
    '

# clean
find groups -mindepth 1 -maxdepth 3 -type d -name "*.raw" | parallel -r rm -rf
find groups -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -rf
find . -mindepth 1 -maxdepth 3 -type f -name "output.*" | parallel -r rm

echo $(find groups -mindepth 1 -maxdepth 1 -type d | wc -l) 
    $(find groups -mindepth 1 -maxdepth 3 -type f -name "*.nwk.pdf" | wc -l)

```

cat taxon/group_target.tsv | sed -e '1d' | grep "^249" |
    parallel --colsep '\t' --no-run-if-empty --line-buffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget:[{4}]\n"
        
        egaz template GENOMES/{2}/{4} \
            $(cat taxon/{2}.sizes | cut -f 1 | grep -v -x "{4}" | xargs -I \[\] echo "GENOMES/{2}/[]") \
            --multi -o groups/{2}/ \
            --order \
            --parallel 24 -v

        bash groups/{2}/1_pair.sh
        bash groups/{2}/3_multi.sh
    '