# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/30 15:57:37
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        return super().parser()
    
    def run(self):
        shell(
            "(pairtools sort"
            " {self.extra}"
            " --nproc {self.snakemake.resources.ncpus}"
            " --memory {self.snakemake.resources.mem_gb}"
            " -o {self.snakemake.output}"
            " {self.snakemake.input}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)