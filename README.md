# KmerGO

KmerGO is a user-friendly tool to identify the group-specific sequences on two groups of high throughput sequencing datasets. A sequence that is present, or rich, in one group, but absent, or scarce, in another group is considered “group-specific” here.

KmerGO offers graphic interface under Windows and Linux with one-click installation free from environmental settings.

Through multi-processes parallel computing, KmerGO is able to capture all the group-specific k-mers (k up to 40bps) for 1.07 TB data in fasta format on a regular standalone server(CPU Intel(R) Xeon(R) E5-2620 v4) in 21.5 hours, including 4 hours Kmer counting with KMC[https://academic.oup.com/bioinformatics/article/33/17/2759/3796399] and 17.5 hours(16 processes number) group-specific kmer identification, and return the assembled group-specific sequences.

## Running on two operation systems

KmerGO can be run directly on Windows and Linux operating systems, without extra enviromental settings or configurations.

### Running on Windows

KmerGO has been tested on Windows 7/8/10 (64 bits).

And you can click [here](https://github.com/ChnMasterOG/KmerGO/releases/download/v1.4.0/KmerGO_for_windows_x64.zip) to download the Windows version of KmerGO.

Decompress the file and run the "Install.bat" with the highest authority if you first use KmerGO.

After that, you can double-click **KmerGO.exe** to run it.

### Running on Linux

KmerGO has been tested on Ubuntu 16, Debian 9, CentOS 7, Fedora 30 and Deepin 15 (64 bits).

And you can click [here](https://github.com/ChnMasterOG/KmerGO/releases/download/v1.4.0/KmerGO_for_linux_x64.zip) to download the Linux version of KmerGO.

Decompress the file and enter the software path.

Type these commands if you first use KmerGO:

> sudo chmod +x KmerGO  
> sudo chmod +x ./bin/*

After that, you can type "**./KmerGO**" to run it.

## Usage

<p align="center">
  <img src="https://raw.githubusercontent.com/ChnMasterOG/KmerGO/master/resource/user_interface.png"/ width="484" height="638">
</p>

**There are four modules in KmerGO**: 1.K-mer counting, 2.K-mer counting union matrix, 3.K-mer features filtering, 4.K-mers assembly. And you can use these four modules to identify the group-specific sequences on two groups or use any module to calculate separately.

### One-Click Running

As shown in the following Gif, You should first select a workpath folder of work path. KmerGO will produce all step results under the work path.

After that, you can set the parameters and press the button "One-Click Start" to start the running. KmerGO can execute all steps automatically if there are no errors.

![alt tag](https://raw.githubusercontent.com/ChnMasterOG/KmerGO/master/resource/oneclick_running.gif)

And the result files will be created in "{$workpath}/contig_result" folder.

### Running modules separately

**1. KMC operation**

If you want obtain the K-mer occurring times using the short reads of FASTA/Q format files independently, you can set it up according to the following steps:

(1) Select "FASTA/Q files path of group A" and "FASTA/Q files path of group B" of input files path.

(2) Select "K-mer counting files path" of output files path.

(3) Set K-mer length (K), K-mer occurring times lower value (MinValue) and K-mer occurring times upper value (MaxValue).

(4) Click "Start" button to run it.

![alt tag](https://raw.githubusercontent.com/ChnMasterOG/KmerGO/master/resource/step_kmc.gif)

**2. Union operation**

If you want obtain the K-mer matrixes in different samples using sorted K-mer occurring times independently, you can set it up according to the following steps:

(1) Select "K-mer counting files path" of input files path.

(2) Select "KMC matrix files path" of output matrixes path.

(3) Set K-mer length (K), K-mer occurring times upper value (MaxValue), the number of samples in two groups and the number of processes.

(4) Click "Start" button to run it.

![alt tag](https://raw.githubusercontent.com/ChnMasterOG/KmerGO/master/resource/step_union.gif)

**3. Filtering operation**

If you want obtain the K-mer filtering features using the K-mer union matrixes independently, you can set it up according to the following steps:

(1) Select "KMC matrix files path" of input matrixes path.

(2) Select "K-mer feature files path" of output features path.

(3) Set the number of samples in two group, ASS threshold value of logical features, rank sum test p value and ASS threshold value of numerical features.

(4) Click "Start" button to run it.

![alt tag](https://raw.githubusercontent.com/ChnMasterOG/KmerGO/master/resource/step_filtering.gif)

**4. CAP3 operation**

If you want obtain sequences assembly using the K-mer features independently, you can set it up according to the following steps:

(1) Select "KMC feature files path" of input features path.

(2) Click "Start" button to run it.

(3) The result files will be created in "contig_result" folder.

![alt tag](https://raw.githubusercontent.com/ChnMasterOG/KmerGO/master/resource/step_cap3.gif)

## Test data

Test data can be found in "release" section. You can also click [here](https://github.com/ChnMasterOG/KmerGO/releases/download/v1.4.0/test_data.zip) to download it.

We prepared two groups of FASTA format files which are stored on "group_A" and "group_B" for END-TO-END testing. 

And we also prepared one testing matrix which is stored on "test_filtering" for running filtering operation independently.

## Contacts and bug reports

Please send bug reports, comments, or questions to

Prof. Ying Wang: [wangying@xmu.edu.cn](mailto:wangying@xmu.edu.cn)

Qi Chen: [23220191151291@stu.xmu.edu.cn](mailto:23220191151291@stu.xmu.edu.cn)

----------

Last update: 9-Jan-2020
