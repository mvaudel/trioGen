### Threading, IO, and job distribution

You can set the number of variants to process in parallel using the command line options. When doing regression analyses, phenotypes are processed in parallel for each variant. If you see that not all cores are used, you can try to include the number of variants processed in parallel, or run different vcf files in parallel. Note, however, that if the performance on your setup is limited by the I/O, this will not have much of an effect. If this is the case please make sure that TrioGen has direct and rapid access to the files (i.e. not accessing them through a network). IO will be dramatically improved by using SSD discs. If all CPUs are running 100%, you can consider running different phenotypes or vcf files on different nodes.

Variant-specific tasks are parallelized when possible using all available resources. You can override the number of threads used by default using the `-Djava.util.concurrent.ForkJoinPool.common.parallelism` argument. 

Example with 32 threads:
```
java -Djava.util.concurrent.ForkJoinPool.common.parallelism=32 -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.Command [parameters]
```

It is possible to run multiple instances of TrioGen from the same _jar_ file, but we recommend using one clean copy of TrioGen per run.


### Memory Settings

The processing might require more than the default memory provided to the Java virtual machine. If you encounter memoyr issues please increase the maximum memory setting using the `Xmx` option. 

Example, for a maximum of 64 GB of memory:
```
java -Xmx64G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.Command [parameters]
```

