using Pkg, PkgTemplates

t = Template(;
           user="ManuelStapper",
           license="MIT",
           authors=["Manuel Stapper <manuel.stapper@googlemail.com>"],
           dir="~/code",
           julia_version=v"0.7",
           ssh=false,
           plugins=[
               TravisCI(),
               Codecov(),
               Coveralls(),
               AppVeyor(),
           ],
       )

generate("Filtering", t)
