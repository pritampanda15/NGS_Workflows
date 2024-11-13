# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 10:37:54
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.java_opts = get_java_opts(self.snakemake)

        self.gvcfs = list(map("--variant {}".format, self.snakemake.input.gvcfs))

        db_action = self.snakemake.params.get("db_action", "create")
        if db_action == "create":
            self.db_action = "--genomicsdb-workspace-path"
        elif db_action == "update":
            self.db_action = "--genomicsdb-update-workspace-path"
        else:
            raise ValueError(
                "invalid option provided to 'params.db_action'; please choose either 'create' or 'update'."
            )

        self.intervals = self.snakemake.input.get("intervals", "")
        if not self.intervals:
            self.intervals = self.snakemake.params.get("intervals", "")

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' GenomicsDBImport"
                " {self.gvcfs}"
                " --intervals {self.intervals}"
                " {self.extra}"
                " --tmp-dir {tmpdir}"
                " {self.db_action} {self.snakemake.output.db}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)