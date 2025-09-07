import hashlib
import os
import sys
import base64
import time
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import pandas as pd
import numpy as np
import subprocess
import tempfile
import re
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import threading
import json
import random
from collections import Counter
import queue


def _obfuscated_license_check():
    _f1 = hashlib.sha256
    _f2 = os.path.exists
    _f3 = "README.md"
    _f4 = "8ed4af899c0eeb135cb548a89e5d2dfa34b0b66fab43851e7b00f215f86a4145"
    _ = sum([ord(c) for c in "license_check"]) % 100

    if not _f2(_f3):
        return False

    try:
        with open(_f3, 'rb') as f:
            _file_data = f.read()
            _actual_hash = _f1(_file_data).hexdigest()
            _expected_hash_parts = [_f4[i:i + 8] for i in range(0, len(_f4), 8)]
            _reconstructed_hash = ''.join(_expected_hash_parts)

            return _actual_hash == _reconstructed_hash
    except:
        return False


def _validate_license():
    try:
        time.sleep(0.1 + random.random() * 0.2)
        return _obfuscated_license_check()
    except:
        return False


def verify_license():
    error_msg = base64.b64decode(
        "UkVBRE1FLm1k5paH5Lu26aqM6K+B5aSx6LSl77yM6K+35qOA5p+l5paH5Lu25piv5ZCm5a2Y5Zyo5LiU5YaF5a655q2j56Gu44CC56iL5bqP5bCG6YCA5Ye6"
    ).decode()

    if not _validate_license():
        try:
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror("License Error", error_msg)
            root.destroy()
        except:
            print(error_msg)

        time.sleep(1.5)
        sys.exit((hash("license_failed") % 100) + 1)

    return True

if not verify_license():
    sys.exit(1)

###########################################################################
# DNA探针分析器UI类 - 用户界面模块
#GitHub：https://github.com/shuomo-os/FISH-analysis-Tool-0C9
#作者：shuomo
#可自行修改代码内容，但确保代码文件所在目录包括对应README.md文件
###########################################################################
class DNAProbeAnalyzerUI:
    def __init__(self, root):
        self.root = root
        self.root.title("FISH探针设计与分析工具 V2.0C9")  # 版本号为V2.0C9
        self.root.geometry("1600x1200")
        self.root.configure(bg='#000000')
        
        # 设置图标（如果有的话）
        try:
            self.root.iconbitmap("dna_icon.ico")
        except:
            pass
            
        self.analyzer = DNAProbeAnalyzer()
        self.designer = RNAProbeDesigner()
        self.pause_flag = threading.Event()
        self.cancel_flag = threading.Event()
        self.ui_update_queue = queue.Queue()
        self.setup_ui()
        self.setup_ui_update_handler()
        
    def setup_ui_update_handler(self):
        """设置UI更新处理器"""
        def check_queue():
            try:
                while True:
                    task, args = self.ui_update_queue.get_nowait()
                    if task == "log_message":
                        self.log_text.insert(tk.END, args[0] + "\n")
                        self.log_text.see(tk.END)
                    elif task == "update_progress":
                        self.progress_var.set(args[0])
                    elif task == "update_design_progress":
                        self.design_progress_var.set(args[0])
                    elif task == "update_status":
                        self.status_var.set(args[0])
                    elif task == "update_results":
                        self.result_vars['total'].set(args[0])
                        self.result_vars['valid'].set(args[1])
                        self.result_vars['avg_tm'].set(args[2])
                        self.result_vars['avg_gc'].set(args[3])
                    elif task == "enable_buttons":
                        self.start_button.config(state='normal')
                        self.pause_button.config(state='disabled')
                        self.pause_button.config(text="暂停")
                        self.cancel_button.config(state='disabled')
                    elif task == "update_design_tree":
                        self.update_design_tree(args[0])
                    elif task == "update_design_stats":
                        self.design_probe_count.config(text=args[0])
                        self.design_avg_gc.config(text=args[1])
                        self.design_avg_tm.config(text=args[2])
                    elif task == "enable_design_buttons":
                        self.export_btn.config(state=tk.NORMAL)
                        self.blast_btn.config(state=tk.NORMAL)
                        self.design_btn.config(state=tk.NORMAL)
            except queue.Empty:
                pass
            finally:
                self.root.after(100, check_queue)  # 每100ms检查一次队列
                
        self.root.after(100, check_queue)
        
    def update_design_tree(self, probes):
        """更新设计结果表格"""
        # 清空表格
        for item in self.design_tree.get_children():
            self.design_tree.delete(item)
        
        # 添加数据到表格
        for probe in probes:
            self.design_tree.insert("", tk.END, values=(
                probe['id'],
                probe['sequence'],
                probe['start'],
                probe['end'],
                f"{probe['gc_content']:.2f}",
                f"{probe['tm']:.2f}",
                f"{probe.get('complexity', 0):.3f}",
                probe.get('specificity', '未检查')
            ))
        
    def setup_ui(self):
        """设置用户界面"""
        # 创建主框架
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # 配置网格权重
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        
        # 创建选项卡
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.grid(row=0, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # 分析选项卡
        self.analysis_tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(self.analysis_tab, text="探针分析")
        self.setup_analysis_tab(self.analysis_tab)
        
        # 设计选项卡
        self.design_tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(self.design_tab, text="探针设计")
        self.setup_design_tab(self.design_tab)
        
        # 反向互补工具选项卡
        self.revcomp_tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(self.revcomp_tab, text="反向互补工具")
        self.setup_revcomp_tab(self.revcomp_tab)
        
        # 配置网格权重
        main_frame.rowconfigure(0, weight=1)
        
    def setup_analysis_tab(self, tab):
        """设置分析选项卡"""
        # 标题
        title_label = ttk.Label(tab, text="FISH探针分析工具",
                               font=("Arial", 16, "bold"), foreground="darkblue")
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))
        
        # 输入文件选择
        ttk.Label(tab, text="输入文件:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.input_file_var = tk.StringVar()
        ttk.Entry(tab, textvariable=self.input_file_var, width=50).grid(row=1, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(tab, text="浏览...", command=self.browse_input_file).grid(row=1, column=2, padx=5)
        
        # Tm计算方法
        ttk.Label(tab, text="Tm计算方法:").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.tm_method_var = tk.StringVar(value="santalucia")
        tm_combo = ttk.Combobox(tab, textvariable=self.tm_method_var,
                               values=["santalucia", "wallace", "gc", "nn"], state="readonly", width=47)
        tm_combo.grid(row=2, column=1, sticky=tk.W, padx=5)
        
        # 输出文件
        ttk.Label(tab, text="输出文件:").grid(row=3, column=0, sticky=tk.W, pady=5)
        self.output_file_var = tk.StringVar(value="probe_analysis_results.csv")
        ttk.Entry(tab, textvariable=self.output_file_var, width=50).grid(row=3, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(tab, text="浏览...", command=self.browse_output_file).grid(row=3, column=2, padx=5)
        
        # 输出格式选项
        output_options_frame = ttk.Frame(tab)
        output_options_frame.grid(row=4, column=0, columnspan=3, sticky=tk.W, pady=5)
        
        self.overwrite_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(output_options_frame, text="覆盖现有文件", 
                       variable=self.overwrite_var).pack(side=tk.LEFT, padx=5)
        
        # BLAST设置框架
        blast_frame = ttk.LabelFrame(tab, text="本地BLAST设置")
        blast_frame.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10, padx=5)
        blast_frame.columnconfigure(1, weight=1)
        
        # BLAST程序路径
        ttk.Label(blast_frame, text="BLAST程序路径:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.blast_path_var = tk.StringVar(value="D:\\app\\blast-2.17.0+\\bin\\blastn.exe")
        ttk.Entry(blast_frame, textvariable=self.blast_path_var, width=50).grid(row=0, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(blast_frame, text="浏览...", command=self.browse_blast_path).grid(row=0, column=2, padx=5)
        
        # BLAST数据库路径
        ttk.Label(blast_frame, text="BLAST数据库路径:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.db_path_var = tk.StringVar(value="E:\\Blast_db\\refseq_rna")
        ttk.Entry(blast_frame, textvariable=self.db_path_var, width=50).grid(row=1, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(blast_frame, text="浏览...", command=self.browse_db_path).grid(row=1, column=2, padx=5)
        
        # BLAST输出文件
        ttk.Label(blast_frame, text="BLAST输出文件:").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.blast_output_var = tk.StringVar(value="blast_results.txt")
        ttk.Entry(blast_frame, textvariable=self.blast_output_var, width=50).grid(row=2, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(blast_frame, text="浏览...", command=self.browse_blast_output).grid(row=2, column=2, padx=5)
        
        # 运行BLAST选项
        self.run_blast_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(blast_frame, text="运行BLAST分析", variable=self.run_blast_var).grid(row=3, column=0, columnspan=3, sticky=tk.W, pady=5)
        
        # 控制按钮
        button_frame = ttk.Frame(tab)
        button_frame.grid(row=6, column=0, columnspan=3, pady=20)
        
        self.start_button = ttk.Button(button_frame, text="开始分析", command=self.start_analysis)
        self.start_button.pack(side=tk.LEFT, padx=5)
        
        self.pause_button = ttk.Button(button_frame, text="暂停", command=self.toggle_pause, state=tk.DISABLED)
        self.pause_button.pack(side=tk.LEFT, padx=5)
        
        self.cancel_button = ttk.Button(button_frame, text="取消", command=self.cancel_analysis, state=tk.DISABLED)
        self.cancel_button.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(button_frame, text="清除日志", command=self.clear_log).pack(side=tk.LEFT, padx=5)
        
        # 进度条
        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(tab, variable=self.progress_var, maximum=100)
        self.progress_bar.grid(row=7, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)
        
        # 状态标签
        self.status_var = tk.StringVar(value="就绪")
        status_label = ttk.Label(tab, textvariable=self.status_var)
        status_label.grid(row=8, column=0, columnspan=3, pady=5)
        
        # 日志输出
        ttk.Label(tab, text="分析日志:").grid(row=9, column=0, sticky=tk.W, pady=5)
        self.log_text = scrolledtext.ScrolledText(tab, width=100, height=15, font=("Consolas", 10))
        self.log_text.grid(row=10, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        
        # 结果统计框架
        result_frame = ttk.LabelFrame(tab, text="分析结果统计")
        result_frame.grid(row=11, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)
        result_frame.columnconfigure(1, weight=1)
        
        self.result_vars = {
            'total': tk.StringVar(value="0"),
            'valid': tk.StringVar(value="0"),
            'avg_tm': tk.StringVar(value="0.0"),
            'avg_gc': tk.StringVar(value="0.0")
        }
        
        ttk.Label(result_frame, text="总序列数:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        ttk.Label(result_frame, textvariable=self.result_vars['total']).grid(row=0, column=1, sticky=tk.W, padx=5)
        
        ttk.Label(result_frame, text="有效序列数:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        ttk.Label(result_frame, textvariable=self.result_vars['valid']).grid(row=1, column=1, sticky=tk.W, padx=5)
        
        ttk.Label(result_frame, text="平均Tm值:").grid(row=2, column=0, sticky=tk.W, padx=5, pady=2)
        ttk.Label(result_frame, textvariable=self.result_vars['avg_tm']).grid(row=2, column=1, sticky=tk.W, padx=5)
        
        ttk.Label(result_frame, text="平均GC含量:").grid(row=3, column=0, sticky=tk.W, padx=5, pady=2)
        ttk.Label(result_frame, textvariable=self.result_vars['avg_gc']).grid(row=3, column=1, sticky=tk.W, padx=5)
        
        # 配置网格权重
        tab.columnconfigure(1, weight=1)
        tab.rowconfigure(10, weight=1)

    def setup_design_tab(self, tab):
        """设置探针设计选项卡"""
        # 标题
        title_label = ttk.Label(tab, text="RNA FISH探针设计工具",
                                font=("Arial", 16, "bold"), foreground="darkgreen")
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))

        # 序列输入组
        seq_group = ttk.LabelFrame(tab, text="目标序列")
        seq_group.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5, padx=5)
        seq_group.columnconfigure(0, weight=1)

        self.seq_input = scrolledtext.ScrolledText(seq_group, width=80, height=8, font=("Consolas", 10))
        self.seq_input.grid(row=0, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)
        self.seq_input.insert(tk.END, "请输入目标RNA序列（ATCG格式）或上传FASTA文件(5→3）...")

        load_btn = ttk.Button(seq_group, text="加载FASTA文件", command=self.load_fasta)
        load_btn.grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)

        # 参数设置组
        params_group = ttk.LabelFrame(tab, text="探针设计参数")
        params_group.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5, padx=5)
        params_group.columnconfigure(1, weight=1)

        # 探针长度
        ttk.Label(params_group, text="探针长度:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.probe_length = tk.Spinbox(params_group, from_=15, to=50, width=10)
        self.probe_length.delete(0, tk.END)
        self.probe_length.insert(0, "20")
        self.probe_length.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        # GC含量范围
        ttk.Label(params_group, text="GC含量范围 (%):").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        gc_frame = ttk.Frame(params_group)
        gc_frame.grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)

        self.min_gc = tk.Spinbox(gc_frame, from_=20, to=80, width=5, format="%.1f", increment=0.1)
        self.min_gc.delete(0, tk.END)
        self.min_gc.insert(0, "40.0")
        self.min_gc.pack(side=tk.LEFT)

        ttk.Label(gc_frame, text="到").pack(side=tk.LEFT, padx=5)

        self.max_gc = tk.Spinbox(gc_frame, from_=20, to=80, width=5, format="%.1f", increment=0.1)
        self.max_gc.delete(0, tk.END)
        self.max_gc.insert(0, "60.0")
        self.max_gc.pack(side=tk.LEFT)

        ttk.Label(gc_frame, text="%").pack(side=tk.LEFT, padx=5)

        # 熔解温度范围
        ttk.Label(params_group, text="熔解温度范围 (°C):").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        tm_frame = ttk.Frame(params_group)
        tm_frame.grid(row=1, column=1, sticky=tk.W, padx=5, pady=5)

        self.min_tm = tk.Spinbox(tm_frame, from_=30, to=90, width=5, format="%.1f", increment=0.1)
        self.min_tm.delete(0, tk.END)
        self.min_tm.insert(0, "55.0")
        self.min_tm.pack(side=tk.LEFT)

        ttk.Label(tm_frame, text="到").pack(side=tk.LEFT, padx=5)

        self.max_tm = tk.Spinbox(tm_frame, from_=30, to=90, width=5, format="%.1f", increment=0.1)
        self.max_tm.delete(0, tk.END)
        self.max_tm.insert(0, "70.0")
        self.max_tm.pack(side=tk.LEFT)

        ttk.Label(tm_frame, text="°C").pack(side=tk.LEFT, padx=5)

        # Tm计算方法
        ttk.Label(params_group, text="Tm计算方法:").grid(row=1, column=2, sticky=tk.W, padx=5, pady=5)
        self.design_tm_method = tk.StringVar(value="santalucia")
        tm_method_combo = ttk.Combobox(params_group, textvariable=self.design_tm_method,
                                       values=["santalucia", "wallace", "gc", "nn"], state="readonly", width=15)
        tm_method_combo.grid(row=1, column=3, sticky=tk.W, padx=5, pady=5)

        # 探针间距
        ttk.Label(params_group, text="探针间距:").grid(row=2, column=0, sticky=tk.W, padx=5, pady=5)
        self.spacing = tk.Spinbox(params_group, from_=1, to=10, width=10)
        self.spacing.delete(0, tk.END)
        self.spacing.insert(0, "3")
        self.spacing.grid(row=2, column=1, sticky=tk.W, padx=5, pady=5)

        # 复杂度筛选
        ttk.Label(params_group, text="最小复杂度:").grid(row=2, column=2, sticky=tk.W, padx=5, pady=5)
        self.min_complexity = tk.Spinbox(params_group, from_=0.0, to=2.0, width=10, format="%.2f", increment=0.05)
        self.min_complexity.delete(0, tk.END)
        self.min_complexity.insert(0, "0.8")
        self.min_complexity.grid(row=2, column=3, sticky=tk.W, padx=5, pady=5)

        # 特异性检查
        self.check_specificity = tk.BooleanVar(value=True)
        ttk.Checkbutton(params_group, text="检查探针特异性", variable=self.check_specificity).grid(
            row=3, column=3, columnspan=2, sticky=tk.W, padx=5, pady=5)

        # 重复序列筛选
        self.filter_repeats = tk.BooleanVar(value=True)
        ttk.Checkbutton(params_group, text="过滤重复序列", variable=self.filter_repeats).grid(
            row=3, column=2, columnspan=2, sticky=tk.W, padx=5, pady=5)

        # 连续相同碱基筛选
        ttk.Label(params_group, text="最大连续相同碱基长度:").grid(row=3, column=0, sticky=tk.W, padx=5, pady=5)
        self.max_homopolymer_length = tk.Spinbox(params_group, from_=1, to=10, width=10)
        self.max_homopolymer_length.delete(0, tk.END)
        self.max_homopolymer_length.insert(0, "3")
        self.max_homopolymer_length.grid(row=3, column=1, sticky=tk.W, padx=5, pady=5)

        # BLAST设置
        blast_group = ttk.LabelFrame(tab, text="BLAST设置")
        blast_group.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5, padx=5)
        blast_group.columnconfigure(1, weight=1)

        # BLAST程序路径
        ttk.Label(blast_group, text="程序路径:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.design_blast_path_var = tk.StringVar(value="D:\\app\\blast-2.17.0+\\bin\\blastn.exe")
        ttk.Entry(blast_group, textvariable=self.design_blast_path_var, width=20).grid(row=0, column=1,
                                                                                       sticky=(tk.W, tk.E), padx=3)
        ttk.Button(blast_group, text="浏览...", command=self.browse_design_blast_path).grid(row=0, column=2, padx=3)

        # BLAST数据库路径
        ttk.Label(blast_group, text="数据库路径:").grid(row=0, column=3, sticky=tk.W, pady=5)
        self.design_db_path_var = tk.StringVar(value="E:\\Blast_db\\refseq_rna")
        ttk.Entry(blast_group, textvariable=self.design_db_path_var, width=20).grid(row=0, column=4,
                                                                                    sticky=(tk.W, tk.E), padx=3)
        ttk.Button(blast_group, text="浏览...", command=self.browse_design_db_path).grid(row=0, column=5, padx=3)

        # 设计按钮和进度条
        button_frame = ttk.Frame(tab)
        button_frame.grid(row=4, column=0, columnspan=3, pady=10)

        self.design_btn = ttk.Button(button_frame, text="设计探针", command=self.start_design)
        self.design_btn.pack(side=tk.LEFT, padx=5)

        self.export_btn = ttk.Button(button_frame, text="导出结果", command=self.export_design_results, state=tk.DISABLED)
        self.export_btn.pack(side=tk.LEFT, padx=5)

        self.blast_btn = ttk.Button(button_frame, text="运行BLAST检查", command=self.run_design_blast, state=tk.DISABLED)
        self.blast_btn.pack(side=tk.LEFT, padx=5)

        self.design_progress_var = tk.DoubleVar()
        self.design_progress_bar = ttk.Progressbar(button_frame, variable=self.design_progress_var, maximum=100)
        self.design_progress_bar.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # 设计结果表格
        results_frame = ttk.LabelFrame(tab, text="设计结果")
        results_frame.grid(row=5, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5, padx=5)
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(0, weight=1)

        # 创建表格
        columns = ("ID", "序列", "起始位置", "结束位置", "GC含量 (%)", "熔解温度 (°C)", "复杂度", "特异性")
        self.design_tree = ttk.Treeview(results_frame, columns=columns, show="headings", height=10)

        for col in columns:
            self.design_tree.heading(col, text=col)
            self.design_tree.column(col, width=100)

        # 添加滚动条
        scrollbar = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.design_tree.yview)
        self.design_tree.configure(yscrollcommand=scrollbar.set)

        self.design_tree.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))

        # 绑定选择事件
        self.design_tree.bind('<<TreeviewSelect>>', self.on_probe_select)

        # 2. 文字详情（独立区域）
        detail_section = ttk.LabelFrame(tab, text="探针详情")
        detail_section.grid(row=6, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5, padx=5)
        detail_section.columnconfigure(0, weight=1)
        self.probe_detail_text = scrolledtext.ScrolledText(
            detail_section, height=5, font=("Consolas", 9))
        self.probe_detail_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # 3. 位置图（独立区域）
        location_section = ttk.LabelFrame(tab, text="探针在mRNA中的位置")
        location_section.grid(row=7, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5, padx=5)
        location_section.columnconfigure(0, weight=1)
        self.mrna_location_canvas = tk.Canvas(
            location_section, width=800, height=50, bg="white")
        self.mrna_location_canvas.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)



        # 统计信息
        stats_frame = ttk.Frame(results_frame)
        stats_frame.grid(row=3, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)

        ttk.Label(stats_frame, text="探针总数:").pack(side=tk.LEFT, padx=5)
        self.design_probe_count = ttk.Label(stats_frame, text="0")
        self.design_probe_count.pack(side=tk.LEFT, padx=5)

        ttk.Label(stats_frame, text="平均GC含量:").pack(side=tk.LEFT, padx=5)
        self.design_avg_gc = ttk.Label(stats_frame, text="0.0")
        self.design_avg_gc.pack(side=tk.LEFT, padx=5)

        ttk.Label(stats_frame, text="平均Tm:").pack(side=tk.LEFT, padx=5)
        self.design_avg_tm = ttk.Label(stats_frame, text="0.0")
        self.design_avg_tm.pack(side=tk.LEFT, padx=5)

        # 配置网格权重
        tab.columnconfigure(0, weight=1)
        tab.rowconfigure(5, weight=1)

    def on_probe_select(self, event):
        """当选择探针时显示详情"""
        selection = self.design_tree.selection()
        if not selection:
            return

        item = selection[0]
        values = self.design_tree.item(item, 'values')

        if not hasattr(self, 'design_results') or not self.design_results:
            return

        # 查找选中的探针
        probe_id = int(values[0])
        probe = next((p for p in self.design_results if p['id'] == probe_id), None)

        if not probe:
            return

        # 获取原始序列
        original_seq = self.seq_input.get(1.0, tk.END).strip().upper()

        # 显示探针详情
        detail_text = f"探针 ID: {probe['id']}\n"
        detail_text += f"DNA探针序列(5'→3'）: {probe['sequence']}\n"
        detail_text += f"目标RNA片段(5'→3'）: {probe.get('rna_fragment', 'N/A')}\n"
        detail_text += f"位置: {probe['start']}-{probe['end']}\n"
        detail_text += f"GC含量: {probe['gc_content']:.2f}%\n"
        detail_text += f"熔解温度: {probe['tm']:.2f}°C\n"
        detail_text += f"复杂度: {probe.get('complexity', 'N/A'):.3f}\n"
        detail_text += f"特异性: {probe.get('specificity', '未检查')}\n\n"

        # 显示与原始序列的互补配对
        if original_seq and len(original_seq) >= probe['end']:
            # 获取目标RNA片段
            target_rna_seq = original_seq[probe['start'] - 1:probe['end']]
            # 获取DNA探针序列并反向（3'到5'）
            dna_probe_seq = probe['sequence'][::-1]

            detail_text += "DNA探针与目标RNA的互补配对:\n"

            # 显示目标RNA序列
            detail_text += f"目标RNA: 5'-{target_rna_seq}-3'\n"

            # 显示DNA探针序列（3'到5'）
            detail_text += f"DNA探针: 3'-{dna_probe_seq}-5'\n\n"

            # 显示配对情况
            pairing = ""
            for rna_base, dna_base in zip(target_rna_seq, dna_probe_seq):
                # DNA-RNA杂交配对规则
                if (rna_base == 'A' and dna_base == 'T') or \
                        (rna_base == 'T' and dna_base == 'A') or \
                        (rna_base == 'C' and dna_base == 'G') or \
                        (rna_base == 'G' and dna_base == 'C'):
                    pairing += "|"  # 完全配对
                else:
                    pairing += "×"  # 不配对

            detail_text += "配对情况:  " + pairing + "\n"
            detail_text += "           " + "".join([" " if p == " " else "↓" for p in pairing]) + "\n"

            # 计算并显示匹配百分比
            match_count = pairing.count('|')
            total_bases = len(pairing)
            match_percent = (match_count / total_bases) * 100 if total_bases > 0 else 0
            detail_text += f"\n匹配百分比: {match_percent:.1f}% ({match_count}/{total_bases} 碱基匹配)\n"

            # 如果有BLAST结果，显示BLAST信息
            if 'blast_hits' in probe and probe['blast_hits']:
                detail_text += f"\nBLAST命中: {len(probe['blast_hits'])} 个\n"
                for i, hit in enumerate(probe['blast_hits'][:3]):  # 显示前3个命中
                    detail_text += f"{i + 1}. {hit['subject']} (相似度: {hit['identity']}%, E值: {hit['evalue']})\n"

        self.probe_detail_text.delete(1.0, tk.END)
        self.probe_detail_text.insert(tk.END, detail_text)

        # 绘制探针在mRNA中的位置
        self.draw_mrna_location_probes(probe, original_seq)

    def draw_mrna_location_probes(self, probe, original_seq):
        """绘制探针在mRNA中的位置"""
        if not hasattr(self, 'mrna_location_canvas'):
            return

        canvas = self.mrna_location_canvas
        canvas.delete("all")  # 清空画布

        # 设置画布参数
        canvas_width = canvas.winfo_width()
        canvas_height = canvas.winfo_height()
        seq_length = len(original_seq)

        # 绘制mRNA序列的水平线
        mRNA_line_y = canvas_height // 2
        canvas.create_line(50, mRNA_line_y, canvas_width - 50, mRNA_line_y, fill="black", width=2)

        # 绘制探针覆盖区域
        probe_start = probe['start']
        probe_end = probe['end']
        probe_length = probe_end - probe_start + 1

        # 将探针位置映射到画布上
        probe_start_x = 50 + (probe_start - 1) * (canvas_width - 100) / seq_length
        probe_end_x = 50 + (probe_end - 1) * (canvas_width - 100) / seq_length

        # 绘制探针覆盖区域的矩形
        canvas.create_rectangle(probe_start_x, mRNA_line_y - 10, probe_end_x, mRNA_line_y + 10, fill="blue",
                                outline="black")

        # 添加文本标签
        canvas.create_text((probe_start_x + probe_end_x) / 2, mRNA_line_y - 20,
                           text=f"探针{probe['id']} ({probe_start}-{probe_end})", fill="black")
    def setup_revcomp_tab(self, tab):
        """设置反向互补工具选项卡"""
        # 标题
        title_label = ttk.Label(tab, text="基因序列批量反向互补工具",
                               font=("Arial", 16, "bold"), foreground="darkred")
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))
        
        # 文件操作组
        file_group = ttk.LabelFrame(tab, text="文件操作")
        file_group.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5, padx=5)
        file_group.columnconfigure(1, weight=1)
        
        # 文件路径选择
        ttk.Label(file_group, text="CSV文件路径:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.revcomp_file_var = tk.StringVar()
        ttk.Entry(file_group, textvariable=self.revcomp_file_var, width=50).grid(row=0, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(file_group, text="浏览...", command=self.browse_revcomp_file).grid(row=0, column=2, padx=5)
        
        # 加载文件按钮
        self.revcomp_load_btn = ttk.Button(file_group, text="加载文件", command=self.load_revcomp_file, state=tk.DISABLED)
        self.revcomp_load_btn.grid(row=1, column=0, columnspan=3, pady=5)
        
        # 数据处理组
        process_group = ttk.LabelFrame(tab, text="数据处理")
        process_group.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5, padx=5)
        process_group.columnconfigure(1, weight=1)
        
        # 列选择
        ttk.Label(process_group, text="选择序列列:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.revcomp_column_combo = ttk.Combobox(process_group, state="readonly", width=47)
        self.revcomp_column_combo.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=5)
        self.revcomp_column_combo.set("请先加载文件")
        self.revcomp_column_combo.config(state="disabled")
        
        # 操作按钮
        button_frame = ttk.Frame(process_group)
        button_frame.grid(row=1, column=0, columnspan=2, pady=10)
        
        self.revcomp_preview_btn = ttk.Button(button_frame, text="预览结果", command=self.preview_revcomp, state=tk.DISABLED)
        self.revcomp_preview_btn.pack(side=tk.LEFT, padx=5)
        
        self.revcomp_process_btn = ttk.Button(button_frame, text="处理并保存", command=self.process_revcomp, state=tk.DISABLED)
        self.revcomp_process_btn.pack(side=tk.LEFT, padx=5)
        
        # 进度条
        self.revcomp_progress_var = tk.DoubleVar()
        self.revcomp_progress_bar = ttk.Progressbar(process_group, variable=self.revcomp_progress_var, maximum=100)
        self.revcomp_progress_bar.grid(row=2, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        self.revcomp_progress_bar.grid_remove()  # 初始隐藏
        
        # 结果显示组
        result_group = ttk.LabelFrame(tab, text="结果预览")
        result_group.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5, padx=5)
        result_group.columnconfigure(0, weight=1)
        result_group.rowconfigure(0, weight=1)
        
        # 创建分割面板
        paned_window = ttk.PanedWindow(result_group, orient=tk.HORIZONTAL)
        paned_window.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)
        
        # 原始数据预览
        original_frame = ttk.Frame(paned_window)
        ttk.Label(original_frame, text="原始数据:").pack(anchor=tk.W)
        self.original_text = scrolledtext.ScrolledText(original_frame, width=40, height=15, font=("Consolas", 9))
        self.original_text.pack(fill=tk.BOTH, expand=True)
        paned_window.add(original_frame, weight=1)
        
        # 处理结果预览
        result_frame = ttk.Frame(paned_window)
        ttk.Label(result_frame, text="处理结果:").pack(anchor=tk.W)
        self.revcomp_result_text = scrolledtext.ScrolledText(result_frame, width=40, height=15, font=("Consolas", 9))
        self.revcomp_result_text.pack(fill=tk.BOTH, expand=True)
        paned_window.add(result_frame, weight=1)
        
        # 配置网格权重
        tab.columnconfigure(0, weight=1)
        tab.rowconfigure(3, weight=1)
        
    def browse_revcomp_file(self):
        """浏览反向互补工具输入文件"""
        filename = filedialog.askopenfilename(
            title="选择CSV文件",
            filetypes=[("CSV文件", "*.csv"), ("所有文件", "*.*")]
        )
        if filename:
            self.revcomp_file_var.set(filename)
            self.revcomp_load_btn.config(state=tk.NORMAL)
            
    def load_revcomp_file(self):
        """加载反向互补工具文件"""
        file_path = self.revcomp_file_var.get()
        if not os.path.exists(file_path):
            messagebox.showerror("错误", "文件不存在！")
            return
            
        try:
            self.revcomp_df = pd.read_csv(file_path)
            self.revcomp_current_file = file_path
            
            # 更新列选择框
            self.revcomp_column_combo.config(state="normal")
            self.revcomp_column_combo['values'] = self.revcomp_df.columns.tolist()
            if self.revcomp_df.columns.tolist():
                self.revcomp_column_combo.current(0)
            self.revcomp_column_combo.config(state="readonly")
            
            # 启用按钮
            self.revcomp_preview_btn.config(state=tk.NORMAL)
            self.revcomp_process_btn.config(state=tk.NORMAL)
            
            # 显示原始数据预览
            self.show_revcomp_data_preview()
            
            self.log_message(f'成功加载文件: {os.path.basename(file_path)}，共 {len(self.revcomp_df)} 行')
            
        except Exception as e:
            messagebox.showerror("错误", f'加载文件时出错:\n{str(e)}')
            
    def show_revcomp_data_preview(self):
        """显示反向互补工具的原始数据预览"""
        if hasattr(self, 'revcomp_df') and self.revcomp_df is not None:
            # 显示前10行数据
            preview_df = self.revcomp_df.head(10)
            preview_text = preview_df.to_string(index=False)
            self.original_text.delete(1.0, tk.END)
            self.original_text.insert(tk.END, preview_text)
            
    def get_reverse_complement(self, dna_sequence):
        """计算反向互补序列"""
        if pd.isna(dna_sequence) or not isinstance(dna_sequence, str):
            return dna_sequence
        try:
            return str(Seq(dna_sequence).reverse_complement())
        except Exception as e:
            return f"错误: {str(e)}"

    def preview_revcomp(self):
        """预览反向互补结果"""
        if not hasattr(self, 'revcomp_df') or self.revcomp_df is None:
            return

        selected_column = self.revcomp_column_combo.get()

        # 创建副本进行处理预览
        preview_df = self.revcomp_df.head(10).copy()
        preview_df['反向互补序列'] = preview_df[selected_column].apply(self.get_reverse_complement)

        # 显示结果
        preview_text = preview_df.to_string(index=False)
        self.revcomp_result_text.delete(1.0, tk.END)
        self.revcomp_result_text.insert(tk.END, preview_text)

        self.log_message('反向互补结果预览已生成')

        
    def process_revcomp(self):
        """处理并保存反向互补结果"""
        if not hasattr(self, 'revcomp_df') or self.revcomp_df is None:
            return
            
        selected_column = self.revcomp_column_combo.get()
        
        # 选择保存文件路径
        save_path = filedialog.asksaveasfilename(
            title='保存结果', 
            defaultextension='.csv',
            filetypes=[('CSV文件', '*.csv')]
        )
        
        if not save_path:
            return
            
        # 显示进度条
        self.revcomp_progress_bar.grid()
        self.revcomp_progress_bar['value'] = 0
        self.revcomp_progress_bar['maximum'] = len(self.revcomp_df)
        
        try:
            # 处理数据
            result_df = self.revcomp_df.copy()
            for i, row in enumerate(self.revcomp_df.iterrows()):
                result_df.at[row[0], '反向互补序列'] = self.get_reverse_complement(row[1][selected_column])
                self.revcomp_progress_bar['value'] = i + 1
                self.root.update_idletasks()  # 更新UI
                
            # 保存结果
            result_df.to_csv(save_path, index=False, encoding='utf-8-sig')
            
            # 隐藏进度条
            self.revcomp_progress_bar.grid_remove()
            
            messagebox.showinfo("成功", f'处理完成！文件已保存至:\n{save_path}')
            self.log_message(f'反向互补处理完成，共处理 {len(self.revcomp_df)} 条序列')
            
        except Exception as e:
            self.revcomp_progress_bar.grid_remove()
            messagebox.showerror("错误", f'处理过程中出错:\n{str(e)}')
            
    def browse_input_file(self):
        """浏览输入文件"""
        filename = filedialog.askopenfilename(
            title="选择输入文件",
            filetypes=[("CSV文件", "*.csv"), ("TSV文件", "*.tsv"), ("文本文件", "*.txt"), ("所有文件", "*.*")]
        )
        if filename:
            self.input_file_var.set(filename)
            # 自动设置输出文件名
            base_name = os.path.splitext(os.path.basename(filename))[0]
            self.output_file_var.set(f"{base_name}_analysis_results.csv")
            
    def browse_output_file(self):
        """浏览输出文件"""
        filename = filedialog.asksaveasfilename(
            title="选择输出文件",
            defaultextension=".csv",
            filetypes=[("CSV文件", "*.csv"), ("TSV文件", "*.tsv"), ("Excel文件", "*.xlsx"), ("所有文件", "*.*")]
        )
        if filename:
            self.output_file_var.set(filename)
            
    def browse_blast_path(self):
        """浏览BLAST程序路径"""
        filename = filedialog.askopenfilename(
            title="选择BLAST程序",
            filetypes=[("可执行文件", "*.exe"), ("所有文件", "*.*")]
        )
        if filename:
            self.blast_path_var.set(filename)
            
    def browse_db_path(self):
        """浏览BLAST数据库路径"""
        directory = filedialog.askdirectory(
            title="选择BLAST数据库目录"
        )
        if directory:
            self.db_path_var.set(directory)
            
    def browse_blast_output(self):
        """浏览BLAST输出文件"""
        filename = filedialog.asksaveasfilename(
            title="选择BLAST输出文件",
            defaultextension=".txt",
            filetypes=[("文本文件", "*.txt"), ("所有文件", "*.*")]
        )
        if filename:
            self.blast_output_var.set(filename)
            
    def browse_design_blast_path(self):
        """浏览设计选项卡中的BLAST程序路径"""
        filename = filedialog.askopenfilename(
            title="选择BLAST程序",
            filetypes=[("可执行文件", "*.exe"), ("所有文件", "*.*")]
        )
        if filename:
            self.design_blast_path_var.set(filename)
            
    def browse_design_db_path(self):
        """浏览设计选项卡中的BLAST数据库路径"""
        directory = filedialog.askdirectory(
            title="选择BLAST数据库目录")
        if directory:
            self.design_db_path_var.set(directory)
            
    def load_fasta(self):
        """加载FASTA文件"""
        filepath = filedialog.askopenfilename(
            title="打开FASTA文件", 
            filetypes=[("FASTA文件", "*.fasta *.fa *.fna"), ("所有文件", "*.*")]
        )
        
        if filepath:
            try:
                # 使用简单的文件读取而不是Bio.SeqIO
                with open(filepath, 'r') as f:
                    content = f.read()
                
                # 提取序列部分（简单处理，跳过标题行）
                lines = content.split('\n')
                sequence_lines = [line for line in lines if not line.startswith('>') and line.strip()]
                sequence = ''.join(sequence_lines)
                
                self.seq_input.delete(1.0, tk.END)
                self.seq_input.insert(tk.END, sequence)
                self.log_message(f"已加载序列来自: {os.path.basename(filepath)}")
            except Exception as e:
                messagebox.showerror("错误", f"读取FASTA文件时出错: {str(e)}")
    
    def log_message(self, message):
        """在日志区域添加消息"""
        timestamp = time.strftime("%H:%M:%S")
        self.ui_update_queue.put(("log_message", [f"[{timestamp}] {message}"]))
        
    def clear_log(self):
        """清除日志"""
        self.log_text.delete(1.0, tk.END)
        
    def update_progress(self, value):
        """更新进度条"""
        self.ui_update_queue.put(("update_progress", [value]))
        
    def update_design_progress(self, value):
        """更新设计进度条"""
        self.ui_update_queue.put(("update_design_progress", [value]))
        
    def update_status(self, message):
        """更新状态信息"""
        self.ui_update_queue.put(("update_status", [message]))
        
    def toggle_pause(self):
        """切换暂停状态"""
        if self.pause_flag.is_set():
            self.pause_flag.clear()
            self.pause_button.config(text="暂停")
            self.log_message("▶️ 分析继续")
        else:
            self.pause_flag.set()
            self.pause_button.config(text="继续")
            self.log_message("⏸️ 分析暂停")
            
    def cancel_analysis(self):
        """取消分析"""
        self.cancel_flag.set()
        self.log_message("⏹️ 分析取消请求已发送")
        
    def start_analysis(self):
        """开始分析"""
        input_file = self.input_file_var.get()
        if not input_file:
            messagebox.showerror("错误", "请选择输入文件")
            return
            
        if not os.path.exists(input_file):
            messagebox.showerror("错误", f"输入文件不存在: {input_file}")
            return
            
        # 检查输出文件是否已存在
        output_file = self.output_file_var.get()
        if os.path.exists(output_file) and not self.overwrite_var.get():
            if not messagebox.askyesno("确认", f"输出文件已存在: {output_file}\n是否覆盖?"):
                return
            
        # 重置标志
        self.pause_flag.clear()
        self.cancel_flag.clear()
        
        # 禁用开始按钮，启用暂停和取消按钮
        self.start_button.config(state='disabled')
        self.pause_button.config(state='normal')
        self.cancel_button.config(state='normal')
        
        # 在新线程中运行分析，避免界面冻结
        thread = threading.Thread(target=self.run_analysis)
        thread.daemon = True
        thread.start()
        
    def run_analysis(self):
        """运行分析的主函数"""
        try:
            self.update_status("正在分析...")
            self.log_message("=== 开始DNA探针分析 ===")
            
            # 获取参数
            config = {
                'input_file': self.input_file_var.get(),
                'output_file': self.output_file_var.get(),
                'tm_method': self.tm_method_var.get(),
                'run_blast': self.run_blast_var.get(),
                'blast_path': self.blast_path_var.get(),
                'db_path': self.db_path_var.get(),
                'blast_output': self.blast_output_var.get()
            }
            
            self.log_message(f"输入文件: {config['input_file']}")
            self.log_message(f"Tm计算方法: {config['tm_method']}")
            
            # 读取文件获取总行数
            try:
                # 尝试不同的分隔符
                try:
                    df = pd.read_csv(config['input_file'])
                except:
                    df = pd.read_csv(config['input_file'], sep='\t')
                total_rows = len(df)
            except Exception as e:
                self.log_message(f"读取文件错误: {e}")
                self.ui_update_queue.put(("enable_buttons", []))
                return
                
            # 执行分析
            success = self.analyzer.analyze_probes(
                config=config,
                pause_flag=self.pause_flag,
                cancel_flag=self.cancel_flag,
                progress_callback=self.update_progress,
                log_callback=self.log_message,
                total_rows=total_rows
            )
            
            if self.cancel_flag.is_set():
                self.update_status("分析已取消")
                self.log_message("❌ 分析已取消")
                messagebox.showinfo("信息", "分析已取消")
            elif success:
                self.update_status("分析完成")
                self.log_message("✅ 分析顺利完成！")
                
                # 更新统计信息
                try:
                    results_df = pd.read_csv(config['output_file'])
                    valid_seqs = results_df['valid_sequence'].sum() if 'valid_sequence' in results_df.columns else 0
                    avg_tm = results_df['tm'].mean() if 'tm' in results_df.columns and not results_df['tm'].isnull().all() else 0
                    avg_gc = results_df['gc_content'].mean() if 'gc_content' in results_df.columns and not results_df['gc_content'].isnull().all() else 0
                    
                    self.ui_update_queue.put(("update_results", [
                        str(total_rows),
                        str(valid_seqs),
                        f"{avg_tm:.2f}°C",
                        f"{avg_gc:.2f}%"
                    ]))
                    
                    messagebox.showinfo("完成", f"分析完成！结果已保存到: {config['output_file']}")
                    
                    # 如果选择了运行BLAST，则执行BLAST分析
                    if config['run_blast']:
                        self.run_blast_analysis(config)
                        
                except Exception as e:
                    self.log_message(f"读取结果文件错误: {e}")
                    messagebox.showinfo("完成", f"分析完成！但统计信息生成失败: {e}")
            else:
                self.update_status("分析失败")
                self.log_message("❌ 分析失败")
                messagebox.showerror("错误", "分析过程中出现错误，请查看日志")
                
        except Exception as e:
            self.log_message(f"分析过程中出现未预期错误: {e}")
            messagebox.showerror("错误", f"分析失败: {e}")
        finally:
            self.ui_update_queue.put(("enable_buttons", []))
            self.update_progress(0)
            
    def run_blast_analysis(self, config):
        """运行BLAST分析"""
        try:
            self.update_status("正在运行BLAST分析...")
            self.log_message("=== 开始BLAST分析 ===")
            
            # 检查BLAST程序是否存在
            if not os.path.exists(config['blast_path']):
                self.log_message(f"❌ BLAST程序不存在: {config['blast_path']}")
                messagebox.showerror("错误", f"BLAST程序不存在: {config['blast_path']}")
                return
                
            # 检查数据库是否存在
            db_name = os.path.basename(config['db_path'])
            db_dir = os.path.dirname(config['db_path'])
            db_files = [f for f in os.listdir(db_dir) if f.startswith(db_name)]
            if not db_files:
                self.log_message(f"❌ BLAST数据库不存在: {config['db_path']}")
                messagebox.showerror("错误", f"BLAST数据库不存在: {config['db_path']}")
                return
                
            # 创建临时FASTA文件用于BLAST查询
            temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
            try:
                # 读取分析结果，只选择有效序列
                results_df = pd.read_csv(config['output_file'])
                valid_seqs = results_df[results_df['valid_sequence'] == True]
                
                if len(valid_seqs) == 0:
                    self.log_message("❌ 没有有效的序列可供BLAST分析")
                    messagebox.showwarning("警告", "没有有效的序列可供BLAST分析")
                    return
                    
                # 写入FASTA文件
                for _, row in valid_seqs.iterrows():
                    seq_id = row['id']
                    sequence = row['sequence']
                    temp_fasta.write(f">{seq_id}\n{sequence}\n")
                temp_fasta.close()
                
                self.log_message(f"创建临时FASTA文件: {temp_fasta.name}")
                self.log_message(f"有效序列数: {len(valid_seqs)}")
                
                # 构建BLAST命令 - 添加-evalue 0.1参数
                blast_cmd = [
                    config['blast_path'],
                    "-db", config['db_path'],
                    "-query", temp_fasta.name,
                    "-outfmt", "7",
                    "-max_target_seqs", "30",  # 限制每个查询的最大命中数为30
                    "-evalue", "0.1",  # 添加最低E值参数
                    "-out", config['blast_output']
                ]
                
                self.log_message(f"运行BLAST命令: {' '.join(blast_cmd)}")
                
                # 运行BLAST
                process = subprocess.Popen(
                    blast_cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True
                )
                
                # 等待进程完成
                stdout, stderr = process.communicate()
                
                if process.returncode == 0:
                    self.log_message("✅ BLAST分析完成")
                    self.log_message(f"结果已保存到: {config['blast_output']}")
                    
                    # 合并BLAST结果到Tm分析结果
                    merged_file = config['output_file'].replace('.csv', '_with_blast.csv')
                    if self.merge_blast_results(config['output_file'], config['blast_output'], merged_file):
                        self.log_message(f"✅ 合并结果已保存到: {merged_file}")
                        messagebox.showinfo("完成", f"BLAST分析完成！合并结果已保存到: {merged_file}")
                    else:
                        messagebox.showinfo("完成", f"BLAST分析完成！但合并结果失败，原始结果已保存到: {config['blast_output']}")
                else:
                    self.log_message(f"❌ BLAST分析失败，错误信息: {stderr}")
                    messagebox.showerror("错误", f"BLAST分析失败: {stderr}")
                    
            finally:
                # 删除临时文件
                os.unlink(temp_fasta.name)
                
        except Exception as e:
            self.log_message(f"BLAST分析过程中出现错误: {e}")
            messagebox.showerror("错误", f"BLAST分析失败: {e}")
            
    def merge_blast_results(self, tm_results_file, blast_results_file, output_file):
        """合并Tm分析结果和BLAST结果（提取前20个BLAST命中，但只保留前5个的详细信息，并统计命中总数）"""
        try:
            # 读取Tm分析结果
            tm_df = pd.read_csv(tm_results_file)
            
            # 解析BLAST结果文件，提取每个查询的命中信息
            blast_data = {}
            current_query = None
            current_hits = []
            hit_counts = {}  # 存储每个查询的命中总数
            
            with open(blast_results_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    # 检测新的查询开始
                    if line.startswith('# Query:'):
                        # 保存前一个查询的结果
                        if current_query and current_hits:
                            blast_data[current_query] = current_hits[:5]  # 只保留前5个命中的详细信息
                            hit_counts[current_query] = len(current_hits)  # 保存命中总数
                        
                        # 开始新的查询
                        current_query = line.split()[2]
                        current_hits = []
                    
                    # 跳过注释行和空行
                    elif line.startswith('#') or not line:
                        continue
                    
                    # 处理数据行
                    else:
                        parts = line.split('\t')
                        if len(parts) >= 12:  # 确保有足够的字段
                            hit_info = {
                                'subject': parts[1],  # subject acc.ver
                                'identity': parts[2],  # % identity
                                'evalue': parts[10],   # evalue
                                'bitscore': parts[11]  # bit score
                            }
                            current_hits.append(hit_info)
                
                # 保存最后一个查询的结果
                if current_query and current_hits:
                    blast_data[current_query] = current_hits[:5]  # 只保留前5个命中的详细信息
                    hit_counts[current_query] = len(current_hits)  # 保存命中总数
            
            # 为每个查询创建BLAST结果列
            # 首先添加命中总数列
            tm_df['blast_hits_count'] = tm_df['id'].map(lambda x: hit_counts.get(x, 0))
            
            # 为前5个命中创建详细信息列
            for i in range(1, 6):  # 为前5个命中创建列
                tm_df[f'blast_hit_{i}'] = ''
                tm_df[f'blast_identity_{i}'] = ''
                tm_df[f'blast_evalue_{i}'] = ''
                tm_df[f'blast_bitscore_{i}'] = ''
            
            # 填充BLAST结果
            for idx, row in tm_df.iterrows():
                query_id = row['id']
                if query_id in blast_data:
                    hits = blast_data[query_id]
                    for i, hit in enumerate(hits):
                        if i < 5:  # 确保不超过5个
                            tm_df.at[idx, f'blast_hit_{i+1}'] = hit['subject']
                            tm_df.at[idx, f'blast_identity_{i+1}'] = hit['identity']
                            tm_df.at[idx, f'blast_evalue_{i+1}'] = hit['evalue']
                            tm_df.at[idx, f'blast_bitscore_{i+1}'] = hit['bitscore']
            
            # 保存合并后的结果
            tm_df.to_csv(output_file, index=False)
            return True
            
        except Exception as e:
            self.log_message(f"合并BLAST结果时出错: {e}")
            import traceback
            self.log_message(traceback.format_exc())
            return False
            
    def start_design(self):
        """开始设计探针"""
        sequence = self.seq_input.get(1.0, tk.END).strip().upper()
        
        # 验证序列
        if not sequence or sequence == "请输入目标RNA序列（ATCG格式）或上传FASTA文件...":
            messagebox.showerror("错误", "请输入目标序列")
            return
        
        valid_chars = set('ATCGU')
        if any(char not in valid_chars for char in sequence):
            messagebox.showerror("错误", "序列包含无效字符。只允许A,T,C,G,U")
            return
        
        # 收集参数
        parameters = {
            'probe_length': int(self.probe_length.get()),
            'min_gc': float(self.min_gc.get()),
            'max_gc': float(self.max_gc.get()),
            'min_tm': float(self.min_tm.get()),
            'max_tm': float(self.max_tm.get()),
            'spacing': int(self.spacing.get()),
            'min_complexity': float(self.min_complexity.get()),
            'check_specificity': self.check_specificity.get(),
            'filter_repeats': self.filter_repeats.get(),
            'tm_method': self.design_tm_method.get(),
            'max_homopolymer_length': int(self.max_homopolymer_length.get())
        }
        
        # 禁用设计按钮
        self.design_btn.config(state=tk.DISABLED)
        self.design_progress_var.set(0)
        
        # 在新线程中运行设计，避免界面冻结
        thread = threading.Thread(target=self.run_design, args=(sequence, parameters))
        thread.daemon = True
        thread.start()
    
    def run_design(self, sequence, parameters):
        """运行探针设计"""
        try:
            self.log_message("=== 开始RNA FISH探针设计 ===")
            self.log_message(f"目标序列长度: {len(sequence)} bp")
            self.log_message(f"探针长度: {parameters['probe_length']} bp")
            self.log_message(f"GC范围: {parameters['min_gc']}%-{parameters['max_gc']}%")
            self.log_message(f"Tm范围: {parameters['min_tm']}°C-{parameters['max_tm']}°C")
            self.log_message(f"Tm计算方法: {parameters['tm_method']}")
            self.log_message(f"最小复杂度: {parameters['min_complexity']}")
            self.log_message(f"过滤重复序列: {'是' if parameters['filter_repeats'] else '否'}")
            
            # 执行设计
            probes = self.designer.design_probes(
                target_sequence=sequence,
                parameters=parameters,
                progress_callback=self.update_design_progress,
                log_callback=self.log_message
            )
            
            # 更新设计结果
            self.design_results = probes
            self.ui_update_queue.put(("update_design_tree", [probes]))
            
            # 更新统计信息
            if probes:
                avg_gc = sum(p['gc_content'] for p in probes) / len(probes)
                avg_tm = sum(p['tm'] for p in probes) / len(probes)
                self.ui_update_queue.put(("update_design_stats", [
                    str(len(probes)),
                    f"{avg_gc:.2f}%",
                    f"{avg_tm:.2f}°C"
                ]))
            
            self.ui_update_queue.put(("enable_design_buttons", []))
            self.log_message(f"✅ 探针设计完成，共设计{len(probes)}个探针")
            
        except Exception as e:
            self.log_message(f"❌ 探针设计过程中出错: {str(e)}")
            messagebox.showerror("错误", f"探针设计过程中出错: {str(e)}")
        finally:
            self.design_btn.config(state=tk.NORMAL)
    
    def run_design_blast(self):
        """对设计的探针运行BLAST分析"""
        if not hasattr(self, 'design_results') or not self.design_results:
            messagebox.showerror("错误", "没有设计结果可供BLAST分析")
            return
        
        # 检查BLAST设置
        blast_path = self.design_blast_path_var.get()
        db_path = self.design_db_path_var.get()
        
        if not os.path.exists(blast_path):
            messagebox.showerror("错误", f"BLAST程序不存在: {blast_path}")
            return
        
        # 创建临时FASTA文件
        temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        try:
            # 写入探针序列到临时文件
            for probe in self.design_results:
                temp_fasta.write(f">{probe['id']}\n{probe['sequence']}\n")
            temp_fasta.close()
            
            # 选择BLAST输出文件
            blast_output = filedialog.asksaveasfilename(
                title="保存BLAST结果",
                defaultextension=".txt",
                filetypes=[("文本文件", "*.txt"), ("所有文件", "*.*")]
            )
            
            if not blast_output:
                return
            
            self.log_message("=== 开始BLAST分析设计探针 ===")
            
            # 构建BLAST命令 - 添加-evalue 0.1参数
            blast_cmd = [
                blast_path,
                "-db", db_path,
                "-query", temp_fasta.name,
                "-outfmt", "7",
                "-max_target_seqs", "30",  # 限制每个查询的最大命中数为30
                "-evalue", "0.1",  # 添加最低E值参数
                "-out", blast_output
            ]
            
            self.log_message(f"运行BLAST命令: {' '.join(blast_cmd)}")
            
            # 运行BLAST
            process = subprocess.Popen(
                blast_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
            
            # 等待进程完成
            stdout, stderr = process.communicate()
            
            if process.returncode == 0:
                self.log_message("✅ BLAST分析完成")
                
                # 解析BLAST结果并更新探针特异性
                self.parse_blast_results(blast_output)
                self.ui_update_queue.put(("update_design_tree", [self.design_results]))
                
                messagebox.showinfo("完成", f"BLAST分析完成！结果已保存到: {blast_output}")
            else:
                self.log_message(f"❌ BLAST分析失败，错误信息: {stderr}")
                messagebox.showerror("错误", f"BLAST分析失败: {stderr}")
                
        finally:
            # 删除临时文件
            os.unlink(temp_fasta.name)
    
    def parse_blast_results(self, blast_output_file):
        """解析BLAST结果并更新探针特异性"""
        try:
            blast_results = {}
            
            with open(blast_output_file, 'r') as f:
                current_query = None
                hits = []
                
                for line in f:
                    line = line.strip()
                    
                    if line.startswith('# Query:'):
                        if current_query and hits:
                            blast_results[current_query] = hits
                        
                        current_query = line.split()[2]
                        hits = []
                    
                    elif line.startswith('#') or not line:
                        continue
                    
                    else:
                        parts = line.split('\t')
                        if len(parts) >= 12:
                            hit_info = {
                                'subject': parts[1],
                                'identity': float(parts[2]),
                                'evalue': float(parts[10]),
                                'bitscore': float(parts[11])
                            }
                            hits.append(hit_info)
                
                # 添加最后一个查询的结果
                if current_query and hits:
                    blast_results[current_query] = hits
            
            # 更新探针特异性
            for probe in self.design_results:
                probe_id = str(probe['id'])
                if probe_id in blast_results:
                    hits = blast_results[probe_id]
                    hit_count = len(hits)
                    
                    # 计算特异性评分（基于最高identity、evalue和命中数量）
                    if hits:
                        best_hit = hits[0]
                        identity = best_hit['identity']
                        evalue = best_hit['evalue']
                        
                        # 新的特异性判断逻辑
                        if identity > 99 and evalue < 1e-5 and hit_count < 5:
                            probe['specificity'] = "高"
                        elif identity > 95 and evalue < 1e-4 and hit_count < 10:
                            probe['specificity'] = "中"
                        elif identity > 75 and evalue < 0.01 and hit_count < 20:
                            probe['specificity'] = "低"
                        else:
                            probe['specificity'] = "不建议使用"
                            
                        # 添加BLAST命中详细信息
                        probe['blast_hits'] = hits[:10]  # 保存前10个命中
                else:
                    probe['specificity'] = "无匹配"
                    
        except Exception as e:
            self.log_message(f"解析BLAST结果时出错: {e}")
    
    def export_design_results(self):
        """导出设计结果"""
        if not hasattr(self, 'design_results') or not self.design_results:
            messagebox.showerror("错误", "没有结果可导出")
            return
        
        filepath = filedialog.asksaveasfilename(
            title="保存设计结果",
            defaultextension=".csv",
            filetypes=[("CSV文件", "*.csv"), ("文本文件", "*.txt")]
        )
        
        if filepath:
            try:
                # 使用pandas导出结果，如果没有安装pandas，则使用csv模块
                try:
                    df = pd.DataFrame(self.design_results)
                    # 处理blast_hits列（如果是列表）
                    if 'blast_hits' in df.columns:
                        df['blast_hits'] = df['blast_hits'].apply(
                            lambda x: '; '.join([f"{h['subject']}({h['identity']}%)" for h in x]) if isinstance(x,
                                                                                                                list) else '')
                    df.to_csv(filepath, index=False)
                except ImportError:
                    # 如果没有pandas，使用csv模块
                    import csv
                    with open(filepath, 'w', newline='') as f:
                        writer = csv.writer(f)
                        # 写入标题
                        headers = ['ID', 'DNA_Sequence', 'RNA_Fragment', 'Start', 'End', 'GC_Content', 'Tm',
                                   'Complexity', 'Specificity']
                        if any('blast_hits' in probe for probe in self.design_results):
                            headers.append('BLAST_Hits')
                        writer.writerow(headers)

                        # 写入数据
                        for probe in self.design_results:
                            row = [
                                probe['id'],
                                probe['sequence'],
                                probe.get('rna_fragment', ''),
                                probe['start'],
                                probe['end'],
                                probe['gc_content'],
                                probe['tm'],
                                probe.get('complexity', 0),
                                probe.get('specificity', 'N/A')
                            ]
                            if 'blast_hits' in probe:
                                hits_str = '; '.join([f"{h['subject']}({h['identity']}%)" for h in probe['blast_hits']])
                                row.append(hits_str)
                            writer.writerow(row)
                
                self.log_message(f"设计结果已导出到: {filepath}")
                messagebox.showinfo("成功", f"设计结果已导出到: {filepath}")
            except Exception as e:
                self.log_message(f"导出结果时出错: {str(e)}")
                messagebox.showerror("错误", f"导出结果时出错: {str(e)}")


###########################################################################
# RNA探针设计器类 - 新增功能模块
###########################################################################
class RNAProbeDesigner:
    def __init__(self):
        self.results = []
    
    def calculate_gc(self, sequence):
        """计算序列的GC含量百分比"""
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

    def calculate_tm(self, sequence, method='santalucia'):
        """使用Bio.SeqUtils.MeltingTemp计算熔解温度"""
        if not sequence or not isinstance(sequence, str):
            return None
            
        sequence = sequence.upper().strip()
        
        # 替换U为T，因为MeltingTemp模块处理的是DNA序列
        sequence = sequence.replace('U', 'T')
        
        if not self.check_sequence_validity(sequence):
            return None
            
        try:
            if method == 'santalucia':
                tm = mt.Tm_NN(Seq(sequence), nn_table=mt.DNA_NN1)
            elif method == 'wallace':
                tm = mt.Tm_Wallace(Seq(sequence))
            elif method == 'gc':
                tm = mt.Tm_GC(Seq(sequence))
            elif method == 'nn':
                tm = mt.Tm_NN(Seq(sequence))
            else:
                tm = mt.Tm_NN(Seq(sequence))
            return round(tm, 2)
        except Exception as e:
            return None
    
    def check_sequence_validity(self, sequence):
        """检查序列有效性"""
        if not sequence:
            return False
            
        sequence = sequence.upper().strip()
        # 替换U为T，因为MeltingTemp模块处理的是DNA序列
        sequence = sequence.replace('U', 'T')
        valid_chars = set('ATCGN')
        return all(char in valid_chars for char in sequence) and len(sequence) > 0
    
    def calculate_complexity(self, sequence):
        """计算序列复杂度（基于序列熵）"""
        if len(sequence) <= 1:
            return 0
            
        # 计算碱基频率
        base_counts = Counter(sequence)
        total = len(sequence)
        
        # 计算熵
        entropy = 0
        for count in base_counts.values():
            p = count / total
            entropy -= p * np.log2(p)
            
        # 标准化到0-2范围（最大熵为2，当4种碱基各占25%时）
        return entropy / 2
    
    def has_repeats(self, sequence, target_sequence, min_repeat_length=6):
        """检查序列是否在目标序列中有重复出现"""
        # 检查序列自身是否有重复
        for i in range(len(sequence) - min_repeat_length + 1):
            substring = sequence[i:i+min_repeat_length]
            if sequence.count(substring) > 1:
                return True
                
        # 检查序列是否在目标序列的其他位置出现（除了它本身的位置）
        occurrences = []
        start_idx = 0
        while True:
            idx = target_sequence.find(sequence, start_idx)
            if idx == -1:
                break
            occurrences.append(idx)
            start_idx = idx + 1
            
        # 如果有多个出现位置，则认为是重复
        return len(occurrences) > 1

    def has_homopolymer(self, sequence, max_homopolymer_length):
        """检查序列中是否存在超过指定长度的连续相同碱基"""
        for base in 'ATCG':
            if base * (max_homopolymer_length + 1) in sequence:
                return True
        return False

    def design_probes(self, target_sequence, parameters, progress_callback=None, log_callback=None):
        """设计RNA FISH探针的核心算法"""
        probes = []
        seq = target_sequence.upper()  # 直接使用字符串
        seq_length = len(seq)

        # 获取参数
        probe_length = parameters['probe_length']
        min_gc = parameters['min_gc']
        max_gc = parameters['max_gc']
        min_tm = parameters['min_tm']
        max_tm = parameters['max_tm']
        spacing = parameters['spacing']
        min_complexity = parameters['min_complexity']
        filter_repeats = parameters['filter_repeats']
        tm_method = parameters.get('tm_method', 'santalucia')
        max_homopolymer_length = parameters.get('max_homopolymer_length', 3)  # 新增参数，默认值为3

        # 沿着序列设计探针
        position = 0
        probe_id = 1
        last_progress = -1

        while position < seq_length - probe_length:
            # 更新进度（每1%更新一次）
            if progress_callback:
                current_progress = int((position / seq_length) * 100)
                if current_progress != last_progress:
                    progress_callback(current_progress)
                    last_progress = current_progress

            # 获取候选探针序列（RNA片段）
            rna_fragment = seq[position:position + probe_length]

            # 计算反向互补序列（DNA探针）
            # 将RNA转换为DNA：U→T，然后计算反向互补
            dna_fragment = rna_fragment.replace('U', 'T')
            candidate_seq = Seq(dna_fragment)
            candidate = str(candidate_seq.reverse_complement())

            # 计算GC含量
            gc_content = self.calculate_gc(candidate)

            # 计算熔解温度
            tm = self.calculate_tm(candidate, method=tm_method)

            # 计算复杂度
            complexity = self.calculate_complexity(candidate)

            # 检查重复序列
            has_repeats = self.has_repeats(rna_fragment, seq) if filter_repeats else False

            # 检查是否存在超过指定长度的连续相同碱基
            has_homopolymer = self.has_homopolymer(candidate, max_homopolymer_length)

            # 检查是否符合条件
            if (tm is not None and gc_content is not None and
                    min_gc <= gc_content <= max_gc and
                    min_tm <= tm <= max_tm and
                    complexity >= min_complexity and
                    not has_repeats and
                    not has_homopolymer):  # 新增条件

                if log_callback:
                    log_callback(
                        f"找到探针 {probe_id}: {candidate} (GC: {gc_content:.2f}%, Tm: {tm:.2f}°C, 复杂度: {complexity:.3f})")

                # 添加到探针列表
                probes.append({
                    'id': probe_id,
                    'sequence': candidate,  # 存储DNA探针序列
                    'rna_fragment': rna_fragment,  # 存储对应的RNA片段
                    'start': position + 1,
                    'end': position + probe_length,
                    'gc_content': gc_content,
                    'tm': tm,
                    'complexity': complexity,
                    'specificity': '未检查'
                })
                probe_id += 1
                position += probe_length + spacing - 1
            else:
                position += 1

        return probes


###########################################################################
# DNA探针分析器类 - 核心功能模块
###########################################################################
class DNAProbeAnalyzer:
    def __init__(self):
        self.results = []
    
    #######################################################################
    # 序列分析模块
    #######################################################################
    def calculate_tm(self, sequence, method='santalucia'):
        """计算DNA序列的Tm值"""
        if not sequence or not isinstance(sequence, str):
            return None
            
        sequence = sequence.upper().strip()
        
        if not self.check_sequence_validity(sequence):
            return None
            
        try:
            if method == 'santalucia':
                tm = mt.Tm_NN(Seq(sequence), nn_table=mt.DNA_NN1)
            elif method == 'wallace':
                tm = mt.Tm_Wallace(Seq(sequence))
            elif method == 'gc':
                tm = mt.Tm_GC(Seq(sequence))
            elif method == 'nn':
                tm = mt.Tm_NN(Seq(sequence))
            else:
                tm = mt.Tm_NN(Seq(sequence))
            return round(tm, 2)
        except Exception as e:
            return None
    
    def calculate_gc_content(self, sequence):
        """计算GC含量"""
        if not sequence or not isinstance(sequence, str):
            return None
            
        sequence = sequence.upper().strip()
        
        if not self.check_sequence_validity(sequence):
            return None
            
        try:
            gc_count = sequence.count('G') + sequence.count('C')
            total_bases = len(sequence)
            
            if total_bases == 0:
                return 0.0
                
            gc_content = (gc_count / total_bases) * 100
            return round(gc_content, 2)
        except Exception as e:
            return None
    
    def check_sequence_validity(self, sequence):
        """检查序列有效性"""
        if not sequence:
            return False
            
        sequence = sequence.upper().strip()
        valid_chars = set('ATCGN')
        return all(char in valid_chars for char in sequence) and len(sequence) > 0
    
    #######################################################################
    # 主分析模块
    #######################################################################
    def analyze_probes(self, config, pause_flag=None, cancel_flag=None, progress_callback=None, 
                      log_callback=None, total_rows=None):
        """主分析函数"""
        input_file = config['input_file']
        if not os.path.exists(input_file):
            if log_callback:
                log_callback(f"错误：输入文件 {input_file} 不存在")
            return False
        
        try:
            # 尝试不同的分隔符
            try:
                df = pd.read_csv(input_file)
            except:
                df = pd.read_csv(input_file, sep='\t')
                
            if log_callback:
                log_callback(f"成功读取文件，共 {len(df)} 条序列")
        except Exception as e:
            if log_callback:
                log_callback(f"读取CSV文件错误: {e}")
            return False
        
        # 检查必要的列
        required_columns = ['sequence']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            if log_callback:
                log_callback(f"CSV文件中缺少必要的列: {missing_columns}")
                log_callback(f"文件列名: {list(df.columns)}")
            return False
        
        results = []
        last_progress = -1
        
        for index, row in df.iterrows():
            # 检查是否取消
            if cancel_flag and cancel_flag.is_set():
                if log_callback:
                    log_callback("分析被用户取消")
                return False
                
            # 检查是否暂停
            if pause_flag and pause_flag.is_set():
                if log_callback:
                    log_callback("分析暂停中...")
                while pause_flag.is_set():
                    time.sleep(0.5)
                    if cancel_flag and cancel_flag.is_set():
                        if log_callback:
                            log_callback("分析被用户取消")
                        return False
            
            sequence = str(row['sequence']).strip()
            probe_id = row.get('id', f"probe_{index+1}")
            
            if log_callback:
                if index % 10 == 0:  # 每10条序列记录一次日志
                    log_callback(f"处理探针 {probe_id} ({index+1}/{len(df)})")
            
            # 更新进度（每1%更新一次）
            if progress_callback and total_rows:
                current_progress = int((index + 1) / total_rows * 100)
                if current_progress != last_progress:
                    progress_callback(current_progress)
                    last_progress = current_progress
            
            # 检查序列有效性
            if not self.check_sequence_validity(sequence):
                if log_callback and index % 10 == 0:  # 减少日志输出
                    log_callback(f"探针 {probe_id}: 序列包含无效字符，跳过")
                result = {
                    'id': probe_id,
                    'sequence': sequence,
                    'valid_sequence': False,
                    'tm': None,
                    'gc_content': None,
                }
                results.append(result)
                continue
            
            # 计算Tm值
            tm = self.calculate_tm(sequence, method=config['tm_method'])
            
            # 计算GC含量
            gc_content = self.calculate_gc_content(sequence)
            
            if log_callback and index % 10 == 0:  # 减少日志输出
                log_callback(f"探针 {probe_id}: Tm={tm}°C, GC={gc_content}%")
            
            result = {
                'id': probe_id,
                'sequence': sequence,
                'valid_sequence': True,
                'tm': tm,
                'gc_content': gc_content,
            }
            
            results.append(result)
        
        # 保存结果
        try:
            results_df = pd.DataFrame(results)
            csv_columns = ['id', 'sequence', 'valid_sequence', 'tm', 'gc_content']
            
            # 确保输出目录存在
            output_dir = os.path.dirname(os.path.abspath(config['output_file']))
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir)
            
            # 根据文件扩展名选择保存格式
            output_ext = os.path.splitext(config['output_file'])[1].lower()
            if output_ext == '.csv':
                results_df[csv_columns].to_csv(config['output_file'], index=False)
            elif output_ext == '.tsv':
                results_df[csv_columns].to_csv(config['output_file'], sep='\t', index=False)
            elif output_ext == '.xlsx':
                results_df[csv_columns].to_excel(config['output_file'], index=False)
            else:
                # 默认保存为CSV
                results_df[csv_columns].to_csv(config['output_file'], index=False)
            
            if log_callback:
                log_callback(f"✅ 分析完成！结果已保存到: {config['output_file']}")
            
            return True
            
        except Exception as e:
            if log_callback:
                log_callback(f"❌ 保存结果错误: {e}")
            return False


###########################################################################
# 主函数模块
###########################################################################
def main():
    # 设置高DPI显示（Windows）
    if os.name == 'nt':
        from ctypes import windll
        windll.shcore.SetProcessDpiAwareness(1)
    
    root = tk.Tk()
    app = DNAProbeAnalyzerUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()