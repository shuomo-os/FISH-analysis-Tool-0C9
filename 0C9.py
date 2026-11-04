import hashlib
import os
import sys
import base64
import time
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext, simpledialog
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
import copy
import shutil
from probe_io import read_probe_table, write_probe_table, read_fasta_records
from blast_support import (
    BlastCancelled, BlastError, MIN_QUERY_LENGTH, build_blast_command,
    write_query_fasta, run_blast, parse_blast_results as parse_blast_output,
    summarize_hits, publish_blast_report,
)


def _obfuscated_license_check():
    _f1 = hashlib.sha256
    _f2 = os.path.exists
    _f3 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "README.md")
    _f4 = "8ed4af899c0eeb135cb548a89e5d2dfa34b0b66fab43851e7b00f215f86a4145"
    _ = sum([ord(c) for c in "license_check"]) % 100

    if not _f2(_f3):
        return False

    try:
        with open(_f3, 'rb') as f:
            # Git may check text files out with CRLF on Windows.
            _file_data = f.read().replace(b'\r\n', b'\n')
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
        self.design_cancel_flag = threading.Event()
        self.ui_update_queue = queue.Queue()
        self.design_results = []
        self.design_source_sequence = ""
        self.analysis_busy = False
        self.design_busy = False
        self._closing = False
        self.worker_threads = []
        self.active_job_paths = {}
        self.setup_ui()
        self.setup_ui_update_handler()
        self.root.protocol("WM_DELETE_WINDOW", self.close_application)
        
    def setup_ui_update_handler(self):
        """Drain worker messages on the Tk thread only."""
        def check_queue():
            if self._closing:
                return
            try:
                while True:
                    task, args = self.ui_update_queue.get_nowait()
                    self._dispatch_ui_task(task, args)
            except queue.Empty:
                pass
            finally:
                if not self._closing:
                    self.root.after(100, check_queue)
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
        ttk.Label(blast_frame, text="BLAST数据库前缀:").grid(row=1, column=0, sticky=tk.W, pady=5)
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
                                       values=["santalucia", "wallace", "gc", "nn", "rna_dna"], state="readonly", width=15)
        tm_method_combo.grid(row=1, column=3, sticky=tk.W, padx=5, pady=5)

        # 探针间距
        ttk.Label(params_group, text="探针间空余碱基数:").grid(row=2, column=0, sticky=tk.W, padx=5, pady=5)
        self.spacing = tk.Spinbox(params_group, from_=0, to=10, width=10)
        self.spacing.delete(0, tk.END)
        self.spacing.insert(0, "3")
        self.spacing.grid(row=2, column=1, sticky=tk.W, padx=5, pady=5)

        # 复杂度筛选
        ttk.Label(params_group, text="最小复杂度:").grid(row=2, column=2, sticky=tk.W, padx=5, pady=5)
        self.min_complexity = tk.Spinbox(params_group, from_=0.0, to=1.0, width=10, format="%.2f", increment=0.05)
        self.min_complexity.delete(0, tk.END)
        self.min_complexity.insert(0, "0.8")
        self.min_complexity.grid(row=2, column=3, sticky=tk.W, padx=5, pady=5)

        ttk.Label(params_group, text="全库匹配需点击“运行BLAST检查”").grid(
            row=4, column=2, columnspan=2, sticky=tk.W, padx=5, pady=5)
        ttk.Label(params_group, text="santalucia/nn: DNA/DNA；rna_dna: RNA/DNA").grid(
            row=4, column=0, columnspan=2, sticky=tk.W, padx=5, pady=5)

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
        ttk.Label(blast_group, text="数据库前缀:").grid(row=0, column=3, sticky=tk.W, pady=5)
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

        self.design_cancel_btn = ttk.Button(button_frame, text="取消", command=self.cancel_design,
                                            state=tk.DISABLED)
        self.design_cancel_btn.pack(side=tk.LEFT, padx=5)

        self.design_progress_var = tk.DoubleVar()
        self.design_progress_bar = ttk.Progressbar(button_frame, variable=self.design_progress_var, maximum=100)
        self.design_progress_bar.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # 设计结果表格
        results_frame = ttk.LabelFrame(tab, text="设计结果")
        results_frame.grid(row=5, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5, padx=5)
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(0, weight=1)

        # 创建表格
        columns = ("ID", "序列", "起始位置", "结束位置", "GC含量 (%)", "熔解温度 (°C)", "复杂度", "BLAST匹配状态")
        self.design_tree = ttk.Treeview(results_frame, columns=columns, show="headings", height=10)

        for col in columns:
            self.design_tree.heading(col, text=col)
            self.design_tree.column(col, width=100)
        self.design_tree.column("序列", width=240)
        self.design_tree.column("BLAST匹配状态", width=220)

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
        original_seq = self.design_source_sequence

        # 显示探针详情
        detail_text = f"探针 ID: {probe['id']}\n"
        detail_text += f"DNA探针序列(5'→3'）: {probe['sequence']}\n"
        detail_text += f"目标RNA片段(5'→3'）: {probe.get('rna_fragment', 'N/A')}\n"
        detail_text += f"位置: {probe['start']}-{probe['end']}\n"
        detail_text += f"GC含量: {probe['gc_content']:.2f}%\n"
        detail_text += f"熔解温度: {probe['tm']:.2f}°C\n"
        detail_text += f"复杂度: {probe.get('complexity', 'N/A'):.3f}\n"
        detail_text += f"Tm方法: {probe.get('tm_method', '未知')}\n"
        detail_text += f"BLAST匹配状态: {probe.get('specificity', '未检查')}\n\n"

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
                        (rna_base in 'TU' and dna_base == 'A') or \
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
                detail_text += f"\nBLAST命中序列数: {probe.get('blast_hits_count', 0)}；比对片段数: {probe.get('blast_hsp_count', 0)}\n"
                detail_text += f"全长完全匹配序列数: {probe.get('blast_full_match_count', 0)}\n"
                if probe.get('blast_hit_limit_reached'):
                    detail_text += "达到返回数量上限，匹配仍可能不完整。\n"
                for i, hit in enumerate(probe['blast_hits'][:3]):  # 显示前3个命中
                    detail_text += f"{i + 1}. {hit['subject']} (相似度: {hit['identity']}%, 查询覆盖率: {hit['query_coverage']}%, E值: {hit['evalue']})\n"

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
        if not seq_length:
            return

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
            self.revcomp_df = read_probe_table(file_path)
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

        try:
            self._assert_paths_available([save_path])
        except ValueError as exc:
            messagebox.showerror("文件正在使用", str(exc), parent=self.root)
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
            write_probe_table(result_df, save_path)
            
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
            filetypes=[("CSV文件", "*.csv"), ("TSV文件", "*.tsv"), ("Excel文件", "*.xlsx"), ("文本文件", "*.txt"), ("所有文件", "*.*")]
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
        self._choose_blast_database(self.db_path_var)
            
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
        self._choose_blast_database(self.design_db_path_var)
            
    def load_fasta(self):
        filepath = filedialog.askopenfilename(
            title="打开FASTA文件",
            filetypes=[("FASTA文件", "*.fasta *.fa *.fna"), ("所有文件", "*.*")])
        if not filepath:
            return
        try:
            records = read_fasta_records(filepath)
            selected = 0
            if len(records) > 1:
                choices = "\n".join(f"{i + 1}. {name}（{len(seq)} nt）"
                                    for i, (name, seq) in enumerate(records[:20]))
                if len(records) > 20:
                    choices += f"\n…共 {len(records)} 条记录"
                number = simpledialog.askinteger(
                    "选择目标序列", f"文件包含多条独立序列，请输入要加载的记录编号：\n{choices}",
                    parent=self.root, minvalue=1, maxvalue=len(records))
                if number is None:
                    return
                selected = number - 1
            name, sequence = records[selected]
            if set(sequence) - set("ACGTU"):
                raise ValueError(f"记录 {name} 含非 A、C、G、T、U 字符，请先确认目标序列。")
            self.seq_input.delete(1.0, tk.END)
            self.seq_input.insert(tk.END, sequence)
            self.log_message(f"已加载 {os.path.basename(filepath)} 的记录 {selected + 1}: {name}（{len(sequence)} nt）")
        except Exception as exc:
            messagebox.showerror("错误", f"读取FASTA文件失败: {exc}", parent=self.root)
    
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
        if self.analysis_busy:
            return
        # Read every Tk variable before starting a worker.
        config = {
            'input_file': self.input_file_var.get().strip(),
            'output_file': self.output_file_var.get().strip(),
            'tm_method': self.tm_method_var.get(),
            'run_blast': self.run_blast_var.get(),
            'blast_path': self.blast_path_var.get().strip(),
            'db_path': self.db_path_var.get().strip(),
            'blast_output': self.blast_output_var.get().strip(),
        }
        try:
            if not os.path.isfile(config['input_file']):
                raise ValueError("请选择存在的输入文件。")
            if not config['output_file']:
                raise ValueError("请选择输出文件。")
            outputs = [config['output_file']]
            if config['run_blast']:
                self._validate_blast_paths(config['blast_path'], config['db_path'])
                if not config['blast_output']:
                    raise ValueError("请选择BLAST输出文件。")
                outputs.extend([config['blast_output'], config['blast_output'] + '.queries.csv',
                                self._merged_output_path(config['output_file'])])
            canonical = [os.path.normcase(os.path.realpath(path))
                         for path in [config['input_file']] + outputs]
            if len(set(canonical)) != len(canonical):
                raise ValueError("输入文件、分析结果、BLAST结果及合并结果必须使用不同路径。")
            if any(os.path.isdir(path) for path in outputs):
                raise ValueError("输出路径必须是文件，不能是目录。")
            existing = [path for path in outputs if os.path.exists(path)]
            if existing and not self.overwrite_var.get():
                if not messagebox.askyesno("确认覆盖", "以下输出文件已存在，是否覆盖？\n" + "\n".join(existing),
                                           parent=self.root):
                    return
            self._reserve_job_paths('analysis', [config['input_file']] + outputs)
        except ValueError as exc:
            messagebox.showerror("参数错误", str(exc), parent=self.root)
            return
        self.pause_flag.clear()
        self.cancel_flag.clear()
        self.analysis_busy = True
        self.progress_var.set(0)
        self.start_button.config(state=tk.DISABLED)
        self.pause_button.config(state=tk.NORMAL)
        self.cancel_button.config(state=tk.NORMAL)
        self._start_worker(self.run_analysis, config)

    def run_analysis(self, config):
        try:
            self.update_status("正在分析...")
            self.log_message("=== 开始DNA探针分析 ===")
            self.log_message(f"输入文件: {config['input_file']}；Tm方法: {config['tm_method']}")
            success = self.analyzer.analyze_probes(
                config, pause_flag=self.pause_flag, cancel_flag=self.cancel_flag,
                progress_callback=self.update_progress, log_callback=self.log_message)
            if self.cancel_flag.is_set():
                raise BlastCancelled("分析已取消。")
            if not success:
                raise ValueError("探针分析失败，请查看日志。")
            results = self.analyzer.results_df.copy(deep=True)
            valid = int((results['valid_sequence'] == True).sum())
            tm = pd.to_numeric(results['tm'], errors='coerce').mean()
            gc = pd.to_numeric(results['gc_content'], errors='coerce').mean()
            self.ui_update_queue.put(("update_results", [str(len(results)), str(valid),
                f"{tm:.2f}°C" if pd.notna(tm) else "—",
                f"{gc:.2f}%" if pd.notna(gc) else "—"]))
            completed_message = f"探针分析结果已保存到: {config['output_file']}"
            if config['run_blast']:
                self.ui_update_queue.put(("analysis_blast_phase", []))
                merged = self.run_blast_analysis(config, results)
                completed_message += f"\nBLAST合并结果: {merged}\n原始BLAST结果: {config['blast_output']}"
            self.update_status("分析完成")
            self.update_progress(100)
            self.log_message("✅ " + completed_message)
            self._show_worker_message("showinfo", "完成", completed_message)
        except BlastCancelled:
            self.update_status("分析已取消")
            self.log_message("分析已取消，BLAST进程已停止；此前已完成的分析文件保留。")
        except BlastError as exc:
            self.update_status("BLAST失败；探针分析结果已保存")
            self.log_message(f"BLAST失败: {exc}")
            self._show_worker_message("showerror", "BLAST失败",
                f"{exc}\n探针分析结果已保存到: {config['output_file']}")
        except Exception as exc:
            self.update_status("分析失败")
            self.log_message(f"分析失败: {exc}")
            self._show_worker_message("showerror", "分析失败", str(exc))
        finally:
            self.ui_update_queue.put(("enable_buttons", []))
            
    def run_blast_analysis(self, config, results_df):
        self.update_status("正在运行BLAST（可取消）...")
        queries = [(position, str(row['sequence']).strip().upper())
                   for position, (_, row) in enumerate(results_df.iterrows())
                   if row['valid_sequence'] == True and len(str(row['sequence']).strip()) >= MIN_QUERY_LENGTH]
        if not queries:
            raise BlastError(f"没有可检索的有效DNA序列（至少 {MIN_QUERY_LENGTH} nt）。")
        skipped = len(results_df) - len(queries)
        if skipped:
            self.log_message(f"有 {skipped} 行无效或长度不足 {MIN_QUERY_LENGTH} nt，标记为未检索。")
        hits, query_map = self._execute_blast(config, results_df, queries, self.cancel_flag)
        if self.cancel_flag.is_set():
            raise BlastCancelled("BLAST已取消。")
        merged = self._merged_output_path(config['output_file'])
        self.merge_blast_results(config['output_file'], config['blast_output'], merged,
                                 query_map=query_map, results_df=results_df, blast_data=hits)
        return merged

    def merge_blast_results(self, tm_results_file, blast_results_file, output_file,
                            query_map=None, results_df=None, blast_data=None):
        if not query_map:
            raise ValueError("BLAST合并缺少内部查询ID与原始行的映射。")
        frame = (read_probe_table(tm_results_file) if results_df is None else results_df.copy(deep=True))
        frame = frame.reset_index(drop=True)
        hits_by_query = (parse_blast_output(blast_results_file, query_map)
                         if blast_data is None else blast_data)
        if set(hits_by_query) != set(query_map):
            raise BlastError("BLAST结果与本次查询列表不一致。")
        if len(set(query_map.values())) != len(query_map):
            raise BlastError("BLAST查询映射包含重复的行位置。")
        for column in list(summarize_hits([])) + ['blast_query_id']:
            frame[column] = pd.Series([''] * len(frame), dtype=object)
        for position, row in frame.iterrows():
            frame.at[position, 'specificity'] = (
                f"未检索（长度不足{MIN_QUERY_LENGTH} nt）" if row.get('valid_sequence') == True
                else "未检索（无效序列）")
        for number in range(1, 6):
            for name in ('hit', 'identity', 'evalue', 'bitscore', 'coverage', 'alignment_length'):
                frame[f'blast_{name}_{number}'] = pd.Series([''] * len(frame), dtype=object)
        for query_id, position in query_map.items():
            if not isinstance(position, int) or not 0 <= position < len(frame):
                raise BlastError("BLAST查询映射的行位置超出分析结果范围。")
            hits = hits_by_query[query_id]
            frame.at[position, 'blast_query_id'] = query_id
            for key, value in summarize_hits(hits).items():
                frame.at[position, key] = value
            # Detail columns show the strongest HSP for each of the first five subjects.
            best_by_subject = {}
            for hit in sorted(hits, key=lambda hit: (-hit['bitscore'], hit['evalue'])):
                best_by_subject.setdefault(hit['subject'], hit)
            for number, hit in enumerate(list(best_by_subject.values())[:5], 1):
                for name, key in [('hit', 'subject'), ('identity', 'identity'), ('evalue', 'evalue'),
                                  ('bitscore', 'bitscore'), ('coverage', 'query_coverage'),
                                  ('alignment_length', 'alignment_length')]:
                    frame.at[position, f'blast_{name}_{number}'] = hit[key]
        write_probe_table(frame, output_file)
        return True
            
    def start_design(self):
        if self.design_busy:
            return
        sequence = ''.join(self.seq_input.get(1.0, tk.END).split()).upper()
        try:
            if not sequence or set(sequence) - set('ACGTU'):
                raise ValueError("请输入只含 A、C、G、T、U 的目标序列（5′→3′）。")
            parameters = self.designer.validate_parameters({
                'probe_length': self.probe_length.get(),
                'min_gc': self.min_gc.get(), 'max_gc': self.max_gc.get(),
                'min_tm': self.min_tm.get(), 'max_tm': self.max_tm.get(),
                'spacing': self.spacing.get(), 'min_complexity': self.min_complexity.get(),
                'filter_repeats': self.filter_repeats.get(),
                'tm_method': self.design_tm_method.get(),
                'max_homopolymer_length': self.max_homopolymer_length.get(),
            })
        except (ValueError, TypeError) as exc:
            messagebox.showerror("参数错误", str(exc), parent=self.root)
            return
        self.design_cancel_flag.clear()
        self.design_busy = True
        self.design_btn.config(state=tk.DISABLED)
        self.blast_btn.config(state=tk.DISABLED)
        self.export_btn.config(state=tk.DISABLED)
        self.design_cancel_btn.config(state=tk.NORMAL)
        self.design_progress_var.set(0)
        self._dispatch_ui_task('design_results', [[], sequence])
        self._start_worker(self.run_design, sequence, parameters)
    
    def run_design(self, sequence, parameters):
        try:
            self.log_message("=== 开始RNA FISH探针设计 ===")
            self.log_message(f"目标长度: {len(sequence)} nt；探针长度: {parameters['probe_length']} nt；"
                             f"Tm方法: {parameters['tm_method']}；探针间空余碱基: {parameters['spacing']}")
            probes = self.designer.design_probes(
                target_sequence=sequence, parameters=parameters,
                progress_callback=self.update_design_progress, log_callback=self.log_message,
                cancel_flag=self.design_cancel_flag)
            if self.design_cancel_flag.is_set():
                self.log_message("探针设计已取消。")
                return
            self.ui_update_queue.put(("design_results", [probes, sequence]))
            self.log_message(f"✅ 探针设计完成，共 {len(probes)} 条；尚未进行全库BLAST检查。")
        except Exception as exc:
            self.log_message(f"探针设计失败: {exc}")
            self._show_worker_message("showerror", "设计失败", str(exc))
        finally:
            self.ui_update_queue.put(("enable_design_buttons", []))
    
    def run_design_blast(self):
        if self.design_busy:
            return
        if not self.design_results:
            messagebox.showerror("错误", "没有设计结果可供BLAST分析。", parent=self.root)
            return
        config = {'blast_path': self.design_blast_path_var.get().strip(),
                  'db_path': self.design_db_path_var.get().strip()}
        try:
            self._validate_blast_paths(config['blast_path'], config['db_path'])
            if any(len(probe['sequence']) < MIN_QUERY_LENGTH for probe in self.design_results):
                raise ValueError(f"当前BLAST设置要求探针长度至少 {MIN_QUERY_LENGTH} nt。")
        except ValueError as exc:
            messagebox.showerror("参数错误", str(exc), parent=self.root)
            return
        output = filedialog.asksaveasfilename(title="保存BLAST结果", defaultextension='.txt',
                                            filetypes=[("文本文件", "*.txt")])
        if not output:
            return
        try:
            if os.path.isdir(output) or os.path.isdir(output + '.queries.csv'):
                raise ValueError("BLAST输出路径必须是文件，不能是目录。")
            self._reserve_job_paths('design_blast', [output, output + '.queries.csv'])
        except ValueError as exc:
            messagebox.showerror("文件正在使用", str(exc), parent=self.root)
            return
        config['blast_output'] = output
        # Clear prior BLAST evidence immediately; a failed rerun must not look current.
        probes = copy.deepcopy(self.design_results)
        for probe in probes:
            for key in [key for key in probe if key.startswith('blast_')]:
                del probe[key]
            probe['blast_hits'] = []
            probe['specificity'] = '检查中'
        self.design_results = copy.deepcopy(probes)
        self.update_design_tree(self.design_results)
        self.probe_detail_text.delete(1.0, tk.END)
        self.design_busy = True
        self.design_cancel_flag.clear()
        self.design_btn.config(state=tk.DISABLED)
        self.blast_btn.config(state=tk.DISABLED)
        self.export_btn.config(state=tk.DISABLED)
        self.design_cancel_btn.config(state=tk.NORMAL)
        self._start_worker(self._run_design_blast, config, probes)

    def parse_blast_results(self, blast_output_file, probes, query_map, blast_data=None):
        hits_by_query = (parse_blast_output(blast_output_file, query_map)
                         if blast_data is None else blast_data)
        if set(hits_by_query) != set(query_map) or set(query_map.values()) != set(range(len(probes))):
            raise BlastError("BLAST结果与本次设计的探针不一致。")
        updated = copy.deepcopy(probes)
        for query_id, position in query_map.items():
            probe = updated[position]
            for key in [key for key in probe if key.startswith('blast_')]:
                del probe[key]
            hits = hits_by_query[query_id]
            probe.update(summarize_hits(hits))
            probe['blast_query_id'] = query_id
            probe['blast_hits'] = hits
        return updated
    
    def export_design_results(self):
        if self.design_busy or not self.design_results:
            return
        filepath = filedialog.asksaveasfilename(title="保存设计结果", defaultextension='.csv',
            filetypes=[("CSV文件", "*.csv"), ("TSV文件", "*.tsv"), ("Excel文件", "*.xlsx")])
        if not filepath:
            return
        try:
            self._assert_paths_available([filepath])
            frame = pd.DataFrame(copy.deepcopy(self.design_results))
            if 'blast_hits' in frame:
                frame['blast_hits'] = frame['blast_hits'].map(
                    lambda hits: json.dumps(hits, ensure_ascii=False) if isinstance(hits, list) else '')
            write_probe_table(frame, filepath)
            self.log_message(f"设计结果已导出到: {filepath}")
            messagebox.showinfo("成功", f"设计结果已导出到: {filepath}", parent=self.root)
        except Exception as exc:
            self.log_message(f"导出失败: {exc}")
            messagebox.showerror("错误", f"导出失败: {exc}", parent=self.root)

    def _dispatch_ui_task(self, task, args):
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
            for key, value in zip(('total', 'valid', 'avg_tm', 'avg_gc'), args):
                self.result_vars[key].set(value)
        elif task == "show_message":
            getattr(messagebox, args[0])(args[1], args[2], parent=self.root)
        elif task == "analysis_blast_phase":
            self.pause_flag.clear()
            self.pause_button.config(state=tk.DISABLED, text="暂停")
        elif task == "enable_buttons":
            self.analysis_busy = False
            self.active_job_paths.pop('analysis', None)
            self.start_button.config(state=tk.NORMAL)
            self.pause_button.config(state=tk.DISABLED, text="暂停")
            self.cancel_button.config(state=tk.DISABLED)
        elif task == "design_results":
            self.design_results, self.design_source_sequence = args
            self.update_design_tree(self.design_results)
            self.probe_detail_text.delete(1.0, tk.END)
            self.mrna_location_canvas.delete("all")
            count = len(self.design_results)
            self.design_probe_count.config(text=str(count))
            self.design_avg_gc.config(text=(
                f"{sum(p['gc_content'] for p in self.design_results) / count:.2f}%" if count else "—"))
            self.design_avg_tm.config(text=(
                f"{sum(p['tm'] for p in self.design_results) / count:.2f}°C" if count else "—"))
        elif task == "design_blast_results":
            self.design_results = args[0]
            self.update_design_tree(self.design_results)
            self.probe_detail_text.delete(1.0, tk.END)
        elif task == "enable_design_buttons":
            self.design_busy = False
            self.active_job_paths.pop('design_blast', None)
            self.design_btn.config(state=tk.NORMAL)
            result_state = tk.NORMAL if self.design_results else tk.DISABLED
            self.export_btn.config(state=result_state)
            self.blast_btn.config(state=result_state)
            self.design_cancel_btn.config(state=tk.DISABLED)


    def _start_worker(self, target, *args):
        self.worker_threads = [thread for thread in self.worker_threads if thread.is_alive()]
        thread = threading.Thread(target=target, args=args, daemon=True)
        self.worker_threads.append(thread)
        thread.start()

    def _assert_paths_available(self, paths):
        canonical = {os.path.normcase(os.path.realpath(path)) for path in paths}
        for reserved in self.active_job_paths.values():
            if canonical & reserved:
                raise ValueError("所选文件正被另一个任务读取或写入，请选择不同路径或等待该任务完成。")
        return canonical

    def _reserve_job_paths(self, job, paths):
        self.active_job_paths[job] = self._assert_paths_available(paths)


    def _show_worker_message(self, kind, title, message):
        self.ui_update_queue.put(("show_message", [kind, title, message]))


    def close_application(self):
        """Cancel work and reap BLAST before closing the Tk event loop."""
        self._closing = True
        self.cancel_flag.set()
        self.design_cancel_flag.set()
        self.pause_flag.clear()
        self.root.withdraw()
        def finish_close():
            if any(thread.is_alive() for thread in self.worker_threads):
                self.root.after(100, finish_close)
            else:
                self.root.destroy()
        finish_close()


    def _choose_blast_database(self, variable):
        path = filedialog.askopenfilename(
            title="选择BLAST核酸数据库索引文件（如 refseq_rna.nal 或 .nin）",
            filetypes=[("BLAST数据库", "*.nal *.nin *.nsq *.nhr *.ndb *.njs"),
                       ("所有文件", "*.*")])
        if path:
            prefix = re.sub(r"(?:\.\d+)?\.(?:nal|nin|nsq|nhr|ndb|njs)$", "", path,
                            flags=re.IGNORECASE)
            variable.set(prefix)


    @staticmethod
    def _validate_blast_paths(blast_path, db_path):
        if not blast_path or not (os.path.isfile(blast_path) or shutil.which(blast_path)):
            raise ValueError(f"BLAST程序不存在: {blast_path}")
        if not db_path.strip():
            raise ValueError("请输入BLAST数据库前缀（不含 .nin 等扩展名）。")
        if os.path.isdir(db_path):
            raise ValueError("数据库路径应包含数据库名称，例如 /data/refseq_rna，不能只选目录。")


    @staticmethod
    def _merged_output_path(output_file):
        stem, extension = os.path.splitext(output_file)
        return stem + "_with_blast" + (extension or ".csv")


    def _execute_blast(self, config, source_frame, queries, cancel_flag):
        output = os.path.abspath(config['blast_output'])
        os.makedirs(os.path.dirname(output), exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.fish_blast_', dir=os.path.dirname(output)) as folder:
            fasta = os.path.join(folder, 'queries.fasta')
            temporary_output = os.path.join(folder, 'results.txt')
            query_map = write_query_fasta(queries, fasta)
            command = build_blast_command(config['blast_path'], config['db_path'], fasta, temporary_output)
            self.log_message("BLAST参数: blastn-short；word_size=7；E-value≤1000；每个查询最多返回1000条序列。")
            self.log_message(f"数据库: {config['db_path']}；本次查询数: {len(query_map)}")
            process = run_blast(command, cancel_flag)
            if process.stderr.strip():
                self.log_message("BLAST提示: " + process.stderr.strip())
            hits = parse_blast_output(temporary_output, query_map)
            mapping = pd.DataFrame([
                {'blast_query_id': query_id, 'source_row': position + 1,
                 'id': source_frame.iloc[position]['id'],
                 'sequence': source_frame.iloc[position]['sequence']}
                for query_id, position in query_map.items()])
            temporary_mapping = os.path.join(folder, 'queries.csv')
            write_probe_table(mapping, temporary_mapping)
            if cancel_flag.is_set():
                raise BlastCancelled("BLAST已取消。")
            # Publish only after the complete report has been parsed successfully.
            publish_blast_report(temporary_output, temporary_mapping, output)
        return hits, query_map


    def cancel_design(self):
        self.design_cancel_flag.set()
        self.log_message("已请求取消探针设计/BLAST检查。")


    def _run_design_blast(self, config, probes):
        try:
            frame = pd.DataFrame(probes)
            queries = [(position, probe['sequence']) for position, probe in enumerate(probes)]
            hits, query_map = self._execute_blast(config, frame, queries, self.design_cancel_flag)
            if self.design_cancel_flag.is_set():
                raise BlastCancelled("BLAST已取消。")
            updated = self.parse_blast_results(config['blast_output'], probes, query_map, hits)
            self.ui_update_queue.put(("design_blast_results", [updated]))
            self.log_message(f"✅ BLAST完成: {config['blast_output']}")
            self._show_worker_message("showinfo", "完成",
                f"BLAST结果: {config['blast_output']}\n查询与原ID对应表: {config['blast_output']}.queries.csv\n"
                "命中序列数和覆盖率已更新，可导出完整设计结果。")
        except BlastCancelled:
            for probe in probes:
                probe['specificity'] = '检查已取消'
            self.ui_update_queue.put(("design_blast_results", [probes]))
            self.log_message("设计探针BLAST检查已取消，检索进程已停止。")
        except Exception as exc:
            for probe in probes:
                probe['specificity'] = '检查失败'
            self.ui_update_queue.put(("design_blast_results", [probes]))
            self.log_message(f"设计探针BLAST失败: {exc}")
            self._show_worker_message("showerror", "BLAST失败", str(exc))
        finally:
            self.ui_update_queue.put(("enable_design_buttons", []))



###########################################################################
# RNA探针设计器类 - 新增功能模块
###########################################################################
class RNAProbeDesigner:
    TM_METHODS = ('santalucia', 'wallace', 'gc', 'nn', 'rna_dna')

    def __init__(self):
        self.results = []
    
    def calculate_gc(self, sequence):
        """仅对确定的 A/T/C/G/U 序列计算 GC 百分比。"""
        if not self.check_sequence_validity(sequence):
            return None
        sequence = sequence.upper().strip().replace('U', 'T')
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100

    def calculate_tm(self, sequence, method='santalucia'):
        """计算 Tm；rna_dna 的 sequence 必须为目标 RNA 的 5′→3′序列。

        其他方法接受 DNA 探针的 5′→3′序列。R_DNA_NN1 使用 RNA 链
        对应的最近邻参数，因此不能把反向互补的 DNA 探针当作 RNA 输入。
        """
        if not isinstance(sequence, str) or method not in self.TM_METHODS:
            return None
        sequence = sequence.upper().strip().replace('U', 'T')
        # 不允许 Biopython 默默移除 N 等未定碱基后计算较短序列的 Tm。
        if len(sequence) < 2 or not self.check_sequence_validity(sequence):
            return None
        try:
            if method == 'santalucia':
                # Allawi & SantaLucia (1997)，而不是 Breslauer 的 DNA_NN1。
                tm = mt.Tm_NN(Seq(sequence), nn_table=mt.DNA_NN3)
            elif method == 'rna_dna':
                tm = mt.Tm_NN(Seq(sequence), nn_table=mt.R_DNA_NN1)
            elif method == 'wallace':
                tm = mt.Tm_Wallace(Seq(sequence))
            elif method == 'gc':
                tm = mt.Tm_GC(Seq(sequence))
            else:
                tm = mt.Tm_NN(Seq(sequence))
            return round(float(tm), 2) if np.isfinite(tm) else None
        except (ValueError, KeyError, TypeError, ZeroDivisionError, OverflowError):
            return None
    
    def check_sequence_validity(self, sequence):
        """检查用于设计和计算的序列是否只含确定的碱基。"""
        if not isinstance(sequence, str):
            return False
        sequence = sequence.upper().strip().replace('U', 'T')
        return bool(sequence) and set(sequence) <= set('ATCG')
    
    def calculate_complexity(self, sequence):
        """按四种碱基的 Shannon 熵计算复杂度，结果范围为 0–1。"""
        if not self.check_sequence_validity(sequence):
            return 0
        sequence = sequence.upper().strip().replace('U', 'T')
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
            
        # 最大熵为 2 bit，因此除以 2 后范围为 0–1。
        return entropy / 2
    
    def has_repeats(self, sequence, target_sequence, min_repeat_length=6):
        """检查重叠的内部重复片段，以及候选片段在目标中的重复位置。"""
        if min_repeat_length < 1:
            raise ValueError('重复片段长度必须至少为 1')
        sequence = sequence.upper().replace('U', 'T')
        target_sequence = target_sequence.upper().replace('U', 'T')
        seen = set()
        for i in range(len(sequence) - min_repeat_length + 1):
            substring = sequence[i:i+min_repeat_length]
            if substring in seen:
                return True
            seen.add(substring)
        if not sequence:
            return False
        first = target_sequence.find(sequence)
        return first >= 0 and target_sequence.find(sequence, first + 1) >= 0

    def has_homopolymer(self, sequence, max_homopolymer_length):
        """检查序列中是否存在超过指定长度的连续相同碱基"""
        sequence = sequence.upper().replace('U', 'T')
        for base in 'ATCG':
            if base * (max_homopolymer_length + 1) in sequence:
                return True
        return False

    def validate_parameters(self, parameters):
        """验证并规范化设计参数，避免无效区间产生误导性的空结果。"""
        result = dict(parameters)
        integer_fields = ('probe_length', 'spacing', 'max_homopolymer_length')
        numeric_fields = integer_fields + ('min_gc', 'max_gc', 'min_tm', 'max_tm', 'min_complexity')
        result.setdefault('max_homopolymer_length', 3)
        for name in numeric_fields:
            try:
                value = float(result[name])
            except (KeyError, TypeError, ValueError):
                raise ValueError(f'参数 {name} 必须是数值') from None
            if isinstance(result[name], (bool, np.bool_)) or not np.isfinite(value):
                raise ValueError(f'参数 {name} 必须是有限数值')
            if name in integer_fields and not value.is_integer():
                raise ValueError(f'参数 {name} 必须是整数')
            result[name] = int(value) if name in integer_fields else value
        if result['probe_length'] < 2:
            raise ValueError('探针长度必须至少为 2 个碱基')
        if result['spacing'] < 0:
            raise ValueError('探针间距（空余碱基数）不能为负数')
        if result['max_homopolymer_length'] < 1:
            raise ValueError('最大连续相同碱基长度必须至少为 1')
        if not 0 <= result['min_gc'] <= result['max_gc'] <= 100:
            raise ValueError('GC 范围必须满足 0 ≤ 最小值 ≤ 最大值 ≤ 100')
        if result['min_tm'] > result['max_tm']:
            raise ValueError('最小 Tm 不能大于最大 Tm')
        if not 0 <= result['min_complexity'] <= 1:
            raise ValueError('最小复杂度必须介于 0 和 1 之间')
        result.setdefault('tm_method', 'santalucia')
        if result['tm_method'] not in self.TM_METHODS:
            raise ValueError(f"不支持的 Tm 方法: {result['tm_method']}")
        result.setdefault('filter_repeats', True)
        return result

    def design_probes(self, target_sequence, parameters, progress_callback=None, log_callback=None,
                      cancel_flag=None):
        """设计RNA FISH探针的核心算法"""
        self.results = []
        probes = []
        if cancel_flag is not None and cancel_flag.is_set():
            return probes
        if not isinstance(target_sequence, str):
            raise ValueError('目标序列必须是文本')
        seq = ''.join(target_sequence.split()).upper()
        if not self.check_sequence_validity(seq):
            raise ValueError('目标序列不能为空，且只能包含 A/T/C/G/U')
        seq_length = len(seq)

        # 获取参数
        parameters = self.validate_parameters(parameters)
        probe_length = parameters['probe_length']
        min_gc = parameters['min_gc']
        max_gc = parameters['max_gc']
        min_tm = parameters['min_tm']
        max_tm = parameters['max_tm']
        spacing = parameters['spacing']
        min_complexity = parameters['min_complexity']
        filter_repeats = parameters['filter_repeats']
        tm_method = parameters.get('tm_method', 'santalucia')
        max_homopolymer_length = parameters['max_homopolymer_length']
        if probe_length > seq_length:
            raise ValueError('目标序列长度不能小于探针长度')

        # 沿着序列设计探针
        position = 0
        probe_id = 1
        last_progress = -1

        while position <= seq_length - probe_length:
            if cancel_flag is not None and cancel_flag.is_set():
                return []
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
            tm_sequence = rna_fragment if tm_method == 'rna_dna' else candidate
            tm = self.calculate_tm(tm_sequence, method=tm_method)

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
                    'tm_method': tm_method,
                    'complexity': complexity,
                    'specificity': '未检查'
                })
                probe_id += 1
                # spacing 表示前一探针末端与后一探针起点之间未覆盖的碱基数。
                position += probe_length + spacing
            else:
                position += 1

        if cancel_flag is not None and cancel_flag.is_set():
            return []
        self.results = probes
        if progress_callback:
            progress_callback(100)
        return probes


###########################################################################
# DNA探针分析器类 - 核心功能模块
###########################################################################
class DNAProbeAnalyzer:
    TM_METHODS = ('santalucia', 'wallace', 'gc', 'nn')

    def __init__(self):
        self.results = []
        self.results_df = None
    
    #######################################################################
    # 序列分析模块
    #######################################################################
    def calculate_tm(self, sequence, method='santalucia'):
        """计算DNA序列的Tm值"""
        if not isinstance(sequence, str) or method not in self.TM_METHODS:
            return None
        sequence = sequence.upper().strip()
        # N 不可交给 Tm_NN 预处理，否则它会计算移除了 N 的较短序列。
        if len(sequence) < 2 or not self.check_sequence_validity(sequence):
            return None
        try:
            if method == 'santalucia':
                tm = mt.Tm_NN(Seq(sequence), nn_table=mt.DNA_NN3)
            elif method == 'wallace':
                tm = mt.Tm_Wallace(Seq(sequence))
            elif method == 'gc':
                tm = mt.Tm_GC(Seq(sequence))
            else:
                tm = mt.Tm_NN(Seq(sequence))
            return round(float(tm), 2) if np.isfinite(tm) else None
        except (ValueError, KeyError, TypeError, ZeroDivisionError, OverflowError):
            return None
    
    def calculate_gc_content(self, sequence):
        """以全部确定的 A/T/C/G 碱基为分母；含未知碱基时不报精确 GC。"""
        if not self.check_sequence_validity(sequence):
            return None
        sequence = sequence.upper().strip()
        gc_count = sequence.count('G') + sequence.count('C')
        return round((gc_count / len(sequence)) * 100, 2)
    
    def check_sequence_validity(self, sequence):
        """是否为非空且没有未定碱基的 DNA 序列。"""
        if not isinstance(sequence, str):
            return False
        sequence = sequence.upper().strip()
        return bool(sequence) and set(sequence) <= set('ATCG')
    
    #######################################################################
    # 主分析模块
    #######################################################################
    def analyze_probes(self, config, pause_flag=None, cancel_flag=None, progress_callback=None, 
                      log_callback=None, total_rows=None):
        """主分析函数"""
        self.results = []
        self.results_df = None
        tm_method = config.get('tm_method', 'santalucia')
        if tm_method not in self.TM_METHODS:
            if log_callback:
                log_callback(f'不支持的 DNA Tm 方法: {tm_method}')
            return False
        input_file = config['input_file']
        if not os.path.exists(input_file):
            if log_callback:
                log_callback(f"错误：输入文件 {input_file} 不存在")
            return False
        
        try:
            df = read_probe_table(input_file)
                
            if log_callback:
                log_callback(f"成功读取文件，共 {len(df)} 条序列")
        except Exception as e:
            if log_callback:
                log_callback(f"读取输入表格错误: {e}")
            return False
        
        # 检查必要的列
        required_columns = ['sequence']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            if log_callback:
                log_callback(f"输入表格中缺少必要的列: {missing_columns}")
                log_callback(f"文件列名: {list(df.columns)}")
            return False
        
        results = []
        last_progress = -1
        
        for index, (_, row) in enumerate(df.iterrows()):
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
                    time.sleep(0.1)
                    if cancel_flag and cancel_flag.is_set():
                        if log_callback:
                            log_callback("分析被用户取消")
                        return False
            
            raw_sequence = row['sequence']
            # 必须先识别缺失值，否则 str(NaN) 会成为看似合法的 "nan"。
            input_sequence = '' if pd.isna(raw_sequence) else str(raw_sequence)
            sequence = input_sequence.strip().upper()
            raw_id = row.get('id', f"probe_{index+1}")
            probe_id = '' if pd.isna(raw_id) else str(raw_id)
            
            if log_callback:
                if index % 10 == 0:  # 每10条序列记录一次日志
                    log_callback(f"处理探针 {probe_id} ({index+1}/{len(df)})")
            
            # 更新进度（每1%更新一次）
            if progress_callback and len(df):
                current_progress = int((index + 1) / len(df) * 100)
                if current_progress != last_progress:
                    progress_callback(current_progress)
                    last_progress = current_progress
            
            alphabet_valid = bool(sequence) and set(sequence) <= set('ATCGN')
            valid_sequence = self.check_sequence_validity(sequence) and len(sequence) >= 2
            if not sequence:
                sequence_error = '序列为空'
            elif not alphabet_valid:
                sequence_error = '序列包含非 DNA 字符；仅接受 A/T/C/G/N'
            elif 'N' in sequence:
                sequence_error = '含未定碱基 N，未计算 Tm 和 GC'
            elif len(sequence) < 2:
                sequence_error = '序列不足 2 个碱基，无法计算 Tm'
            else:
                sequence_error = ''

            tm = self.calculate_tm(sequence, method=tm_method) if valid_sequence else None
            gc_content = self.calculate_gc_content(sequence) if valid_sequence else None
            if valid_sequence and tm is None:
                sequence_error = '所选方法无法计算此序列的 Tm'
            
            if log_callback and index % 10 == 0:  # 减少日志输出
                if sequence_error:
                    log_callback(f"探针 {probe_id}: {sequence_error}")
                else:
                    log_callback(f"探针 {probe_id}: Tm={tm}°C, GC={gc_content}%")
            
            result = {
                'id': probe_id,
                'sequence': sequence,
                'valid_sequence': valid_sequence,
                'tm': tm,
                'gc_content': gc_content,
                'alphabet_valid': alphabet_valid,
                'tm_calculable': tm is not None,
                'tm_method': tm_method,
                'sequence_error': sequence_error,
                'input_sequence': input_sequence,
            }
            
            results.append(result)
        
        # 保存结果
        if cancel_flag is not None and cancel_flag.is_set():
            if log_callback:
                log_callback('分析被用户取消')
            return False
        try:
            columns = ['id', 'sequence', 'valid_sequence', 'tm', 'gc_content',
                       'alphabet_valid', 'tm_calculable', 'tm_method', 'sequence_error',
                       'input_sequence']
            results_df = pd.DataFrame(results, columns=columns)
            write_probe_table(results_df, config['output_file'])
            self.results = results
            self.results_df = results_df
            
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
    verify_license()
    # 设置高DPI显示（Windows）
    if os.name == 'nt':
        from ctypes import windll
        windll.shcore.SetProcessDpiAwareness(1)
    
    root = tk.Tk()
    app = DNAProbeAnalyzerUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
