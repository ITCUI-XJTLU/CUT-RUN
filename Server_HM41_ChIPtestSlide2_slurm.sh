#!/bin/bash
#SBATCH --job-name=hm41_ch2            
#SBATCH --output=./logs/hm41_ch2.log      
#SBATCH --error=./logs/hm41_ch2.log        
#SBATCH --partition=gqi-32c256g,all-12c128g
#SBATCH --ntasks=1                        # 只运行1个任务
#SBATCH --cpus-per-task=5                 # 每个任务使用的CPU核心数
#SBATCH --time=2-00:00:00                   # 最长运行时间
#SBATCH --mem=40G                         # 内存需求

# 加载conda配置
source /home/users/tfcui23/Stat_Gene/GATK/miniconda/etc/profile.d/conda.sh

# 激活你的R环境
conda activate daesc_r

# 输出R版本（可选调试）
Rscript --version

# 运行你的R脚本
Rscript Server_HM41_ChIPtestSlide2.R
