readme.md中常用的命令   
- `rsync`（remote synchronization）远程同步，同步文件/目录，可只传输文件中发生变化的部分   
选项：`-a`（--archive）归档模式进行同步（递归并保留符号连接、文件权限、时间戳、属组信息、属主信息、设备文件和特殊文件）   
    `-v`（--verbose）详细输出    
    `-P`（--partial --progress），保留部分传输的文件（中断后可继续）、显示传输进度   
`rsync -avP ftp.ncbi.nlm.nih.gov::refseq/release/plasmid/ RefSeq/`    
从服务器`ftp.ncbi.nlm.nih.gov`的`refseq/release/plasmid/`目录同步文件到本地的`RefSeq/`目录。`::`是rsync`协议的标识（使用rsync协议传输）   

- `gzip`压缩、解压缩
选项：`-d`decompress解压
    `-c`stdout解压后内容输出到标准输出。结合`>`指定输出保存文件
    `-f`force强制，输出文件已存在则覆盖而不提示

- MinHash含义
用哈希函数将元素映射成数字，从哈希值里找最小值，并只保留最小值，即MinHash签名，重复多次得到多个签名，再比较计算签名的相似度。   
质粒中：将质粒序列拆分为多个k-mer形成一个k-mer集合，再用多个哈希函数处理每个质粒的k-mer集合，生成MinHash签名，比较计算质粒间的相似度   
过滤后保留不相似的质粒（非冗余质粒）

- `tsv-filter`处理tsv（tab-separated values），属于tsv-utils一部分。对tsv文件中数据基于字段值进行过滤
选项：`--le`less than or equal to小于等于
`2:2000`，2表示列，2000为阈值（字段值满足的小于等于的条件）