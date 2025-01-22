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
