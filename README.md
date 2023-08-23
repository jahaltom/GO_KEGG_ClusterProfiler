ClusterProfiler.r: needs a "list" of folders with DESeq2 results and uses script.sh. Uses GeneAnn.txt.

clusterProfiler_Parallele.r: needs a "list" of folders with DESeq2 results and uses script.sh. Uses GeneAnn.txt. To implement use bellow script:

```
cat list | while read i;do
	cat clusterProfiler_Parallele.r | sed "s/FOLDER/$i/g" > clusterProfiler_Parallele.$i.r
	cat script.sh | sed "s/ClusterProfiler.r/clusterProfiler_Parallele.$i.r/g" > script.$i.sh
	sbatch script.$i.sh
done

rm *Group*
```
