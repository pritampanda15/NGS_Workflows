# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/08 19:47:17
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
            "(hicConvertFormat"
            " {self.extra}"
            " --inputFormat {self.snakemake.params.input_format}"
            " --matrices {self.snakemake.input}"
            " --outputFormat {self.snakemake.params.output_format}"
            " --outFileName {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)