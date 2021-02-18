# ingenuous
Based on the algorithms invented by Richard Tjörnhammar (me) in the impetuous-gfa python package. This Julia package is an ongoing effort.

To use a Julia package you will need a proper environment to develop in. I have choosen to work in the Free GNU/Linux system Guix and use an environment defined here: ```https://github.com/richardtjornhammar/scheming-guix/blob/main/manifests/generic_julia.scm```

I then enter that environment :
```
$ guix enviroment --manifest=generic_julia.scm
```
I start Julia :
```
$ julia
```
and enter the Julia Pkg manager by pressing down ']'
```
julia> ]
pkg> add https://github.com/richardtjornhammar/Ingenuous.jl
```
