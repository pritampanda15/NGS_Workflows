# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/01 14:45:13
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
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)
    
    def run(self):
        shell(
            "(meryl histogram"
            " threads={self.snakemake.threads}"
            " {self.extra}"
            " {self.snakemake.input}"
            " > {self.snakemake.output}"
            ") {self.log}"
        )


if __name__ == "__main__":
    wrapper = Wrapper(snakemake)
    wrapper.run()