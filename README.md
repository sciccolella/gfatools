## Getting Started
Compiling of the `example-usage` file. To use the extensions to the struct we need to also include `gfa-ext` files.
```bash
gcc -o eu example-usage.c gfa-io.c gfa-base.c gfa-ed.c kalloc.c gfa-ext.c -lz
```


## Misc
[rGFA]: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
