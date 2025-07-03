# Simulation test for comparing different differential binding analysis softwares

## Requirement
In this docment, we show how to simulate differential binding cases (Histone modification and Transctription factor) to compare different methods.

We have make every simulation process one single function, we also make one differential binding analysis method one single function. All function is placed in the R file ChIPtest_source.R . If you want to use the functions, please make sure you can import all required packages in the R script.


## Key functions 
### Simulate histone modification 
By default, the function will simulate 20,000 sites across the whole genome, where there are 1,000 true sites and 19,000 false sites.

```
simulate_hm(true.width, base.mu, data_path, bins, count_length) 
```
- true.width: The parameter controls the widths of peaks. The default number is 500. 
- base.mu: The parameter controls the heights of peaks. 
- data_path: the path for saving data
- bins: the parameter is used for visualization. If we want to plot the peaks, we need to get count sequences over the peaks. The parameter is used to determined how to calculate each counts. The default value is 15. 
- count_length: the parameter is also used for visualization. It controls the length of count sequences for each site.

### Simulate transcription factor
By default, the function will simulate 20,000 sites across the whole genome, where there are 1,000 true sites and 19,000 false sites.

```
simulate_tf(data_path, fraglen, 
            base.mu, up.mu, down.mu.1, down.mu.2,
            bins_vis, count_length_vis) 
```
- data_path: The path for saving simulated data
- fraglen: The parameter controls the width of peaks in the TF case. The default value is 10 which will simulate very narrow peaks. 
- base.mu: the parameter controls the heights of peaks in the non-differential sites (false sites). 
- up.mu, down.mu.1, down.mu.2: the two parameters control the heights of the peaks at the differential sites. down.mu.1, down.mu.2 should be the same. If the difference between up.mu and down.mu is very large, then the task would become easy. 
- bins_vis: The parameter is used for the visulization. After we simultate the data, we will randomly select 10 sites and drow them out. We need to get the count number over some interval. The parameter controls the length of interval for one count number. The default number is 15. 
- count_length_vis: The parameter is also used for visulization of simulated data. It determines the length of count sequences. The defaul number is 120. 

### Fucntion for scaw
By this function, we can use the software, CSAW, to analyze the simulated data. 

```
run_csaw_differential_analysis(input_data_dir, output_dir, is.tf, window_width, peak_width, iteration_number
) 
```
- input_data_dir: The path of simulated data
- output_dir: The path for saving result files
- is.tf: If we analyze the HM data, it is false; otherwise, it is true. 
- window_width: The default is 150 for histone modification, 10 for histone factor. 
- peak_width: The length of simulated peaks at each sites. The default values are 1,000 for HM and 500 for TF. 
- iteration_number: We may repeat the simulate test for many times. The parameter record the time of simulation test. 

### Function for ChIPtest
By the function, we can use the ChIPtest with the sliding windows to analyze the simulated data. 

```
run_chiptest_slide2_analysis(
    input_data_dir, output_parent_dir, is.tf,  peak_width,
    chiptest_params = list(
      bins_tf = 25, counts_tf = 120,
      filter_tf = 0, window_width_csaw_tf = 150,
      bins_hm = 15, counts_hm = 120,
      filter_hm = 0, window_width_csaw_hm = 150
    )
)
```
- input_data_dir: The path for simulated data
- output_parent_dir: The path for the result files
- is.tf: true for transcription factor; false for histone modification
- peak_width: The length of simulated peaks at each sites. The default values are 1,000 for HM and 500 for TF. 
- chiptest_params: 
    - bins_tf: ChIPtest needs count sequences. The parameter controls the bin size. The parameter is only useful for the TF case. The default value is 25. 
    - counts_tf: ChIPtest needs count seqeunces. The parameter controls the length of count sequences. The parameter is only useful for the TF case. The default value is 120. 
    - filter_tf, filter_hm: the two parameters are 0 by default. you do not need to change it. 
    - window_width_csaw_tf: This pipeline needs CSAW to detect potential sites at first. The parameter controls the window size of the used CSAW. The parameter is only useful for the TF case.
    - bins_hm: ChIPtest needs count sequences. The parameter controls the bin size. The parameter is only useful for the HM case. The default value is 15.
    - counts_hm: ChIPtest needs count seqeunces. The parameter controls the length of count sequences. The parameter is only useful for the HM case. The default value is 120. 
    - window_width_csaw_hm: This pipeline needs CSAW to detect potential sites at first. The parameter controls the window size of the used CSAW. The parameter is only useful for the HM case.

### Function for MACS + ChIPtest
By this method, we first use MACS to detect the potential peaks (make sure you can use the MACS), then we use ChIPtest to do the hypothesis test.

run_macs3_chiptest_pipeline(data_base_dir, output_base_dir,run_tag, count_bin=15, peak_width=500, is_tf = FALSE,source_file_path)

- data_base_dir: The path of simulated data
- output_base_dir: The path for saving result files
- run_tag: record the number of iteration. Because we need to repeat the experiments for many times.
- count_bin: ChIPtest need the count sequences. The parameter control the bin size. 
- peak_width: The width of each peaks. By default, it is 1000 for HM and 500 for TF. 
- is_tf: FALSE for transcription factor; TRUE for histone modification
- source_file_path: The path of the ChIPtest_source.R file. For example, "/path/to/ChIPtest_source.R"

### Function for MACS + Diffbind

By using this function, we first use MACS to detect the potential sites, and then use Diffbind to do the hypothesis test.

```
run_macs3_diffbind_pipeline( 
    data_directory, 
    output_directory_base, 
    iteration_tag, 
    count_bin=15, 
    is_tf = FALSE, 
    peak_width, 
    source_file_path = "/path/to/ChIPtest_source.R")
```
- data_directory: the path for the simulated data
- output_directory_base: the path for the results
- iteration_tag: The number of iteration
- count_bin: The parameter controls the bin size. 
- is_tf: FALSE for the TF case; TRUE for the HM case
- peak_width: the width of peaks at each sites. The default value is 1,000 for the HM and 500 for the TF. 
- source_file_path: The path of the source file. 


### Function for PePr
By this method, we use PePr to analyze the simulated data

run_pepr_pipeline(
    input_data_dir,
    output_dir_base,
    is_tf,
    iteration_tag = "run1"
) 

- input_data_dir: the path of simulated data
- output_dir_base: the path of the final result files
- is_tf: FALSE for the TF case; TRUE for the HM case
- iteration_tag: iteration number

## Completed workflow
In the second section, you can implement each pipeline as one function. It is very convinent. However, if you want to run the whole expriment, you need to simulate the data first and then use the methods to analyze the simulated data. For example, in the file Server_HM1_csaw1.R, we show the codes to perform the whole work. Server_HM1_csaw1_slurm.sh contains the codes for running codes on SLURM. You can find the codes for other methods as well, such as: Server_HM1_ChIPtestSlide2_slurm.sh, Server_HM1_ChIPtestSlide2.sh, Server_HM1_MACS_Diffbind_slurm.sh, Server_HM1_MACS_Diffbind.sh ...... 

HM1, HM4, HM6 ..... are the different simulation setting, for example, the height or the width of the peaks may be different. 











