## Troubleshooting

When encountering an error, please try and troubleshoot the most likely cause of the problem. If the problem persists, please [open an issue](https://github.com/mvaudel/trioGen/issues).


### Commonly encountered issues

#### Could not find or load main class

* Error Message:
```
Error: Could not find or load main class no.uib.triogen.cmd.*
```

* Troubleshooting:
  * Check that the path to the jar file is correct.
  * Check that the version number of the jar file is correct.


#### Insufficient Memory

* Error Message:
```
OpenJDK 64-Bit Server VM warning: INFO: os::commit_memory(0x00007f292a200000, 3997171712, 0) failed; error='Cannot allocate memory' (errno=12)
#
# There is insufficient memory for the Java Runtime Environment to continue.
# Native memory allocation (mmap) failed to map 3997171712 bytes for committing reserved memory.
# An error report file with more information is saved as:
# folder/hs_err_pid5115.log

```
* Troubleshooting:
  * Check that Java **64-Bit** is installed, **not** 32-Bit
  * Check that the memory given to the tool is available on the machine. If multiple processes are running in parallel, please make sure that the sum of memory allocated to each does not exceed the memory available on the machine. 
  * Java might need more memory than what is set via `-Xmx`, see [here](https://stackoverflow.com/questions/15282178/java-using-up-far-more-memory-than-allocated-with-xmx) for details. Try to reduce the amount of memory allocated via `-Xmx` to leave more memory for Java.
  * Try to increase the amount of memory available for the tool.
  
#### An illegal reflective access operation has occurred

* Error Message:
```
WARNING: An illegal reflective access operation has occurred
WARNING: Illegal reflective access by io.airlift.compress.zstd.UnsafeUtil (file:/mnt/work/marc/moba/trioGen/bin/triogen-0.5.0-beta/lib/aircompressor-0.18.jar) to field java.nio.Buffer.address
WARNING: Please consider reporting this to the maintainers of io.airlift.compress.zstd.UnsafeUtil
WARNING: Use --illegal-access=warn to enable warnings of further illegal reflective access operations
WARNING: All illegal access operations will be denied in a future release

```
* Troubleshooting:
  * This warning can be ignored ([source](https://stackoverflow.com/questions/46230413/jdk9-an-illegal-reflective-access-operation-has-occurred-org-python-core-pysys)).
  
#### No StackTrace

* Error Message:
```
java.lang.NullPointerException

```
* Troubleshooting:
  * Add `-XX:-OmitStackTraceInFastThrow` to your command ([Source](https://stackoverflow.com/questions/4659151/recurring-exception-without-a-stack-trace-how-to-reset)).
  
  
#### Airlift _IncompatibleJvmException_

* Error Message:
```
Exception in thread "ImportThread" java.lang.ExceptionInInitializerError
	at io.airlift.compress.zstd.ZstdFrameCompressor.writeMagic(ZstdFrameCompressor.java:57)
	at io.airlift.compress.zstd.ZstdFrameCompressor.compress(ZstdFrameCompressor.java:143)
	at io.airlift.compress.zstd.ZstdCompressor.compress(ZstdCompressor.java:45)
	...
        Suppressed: java.lang.NoClassDefFoundError: Could not initialize class io.airlift.compress.zstd.UnsafeUtil
		at io.airlift.compress.zstd.ZstdFrameCompressor.writeMagic(ZstdFrameCompressor.java:57)
		at io.airlift.compress.zstd.ZstdFrameCompressor.compress(ZstdFrameCompressor.java:143)
		at io.airlift.compress.zstd.ZstdCompressor.compress(ZstdCompressor.java:45)
		...
Caused by: io.airlift.compress.IncompatibleJvmException: Zstandard requires access to java.nio.Buffer raw address field
	at io.airlift.compress.zstd.UnsafeUtil.<clinit>(UnsafeUtil.java:53)
	...

```
* Troubleshooting:
  * Add `--add-opens java.base/java.nio=ALL-UNNAMED` to your command ([Source](https://github.com/airlift/aircompressor/issues/125)).

 
