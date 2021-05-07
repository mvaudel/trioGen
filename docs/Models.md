## Trio Models

TrioGen proposes a selection of {child, mother, father} trio models for trio and haplotype analyses. Please note that it is currently not possible to implement custom models. If this is something important for your research, please [open an issue](https://github.com/mvaudel/trioGen/issues).


### Haplotypes

In the haplotype `h` model, the regression is conducted against the number of transmitted and non-transmitted tested alleles, h, as done by [Chen _et al._](https://doi.org/10.1101/737106).

```
y = �mnt hmnt + �mt hmt + �ft hft + �fnt hfnt + e                                                     (h)
```

Where y represents the phenotypes, _hmnt_ the number of maternal tested alleles non-transmitted to the child, _hmt_ the number of maternal tested alleles transmitted to the child, _hft_ the number of paternal tested alleles transmitted to the child, and _hfnt_ the number of paternal tested alleles non-transmitted to the child.


### Child, mother, and father trio

In the {child, mother, father} `cmf` model, the regression is conducted against the number of tested alleles in the child, mother, and father conditioned against each other.

```
y = �m m + �c c + �f f + e                                                                            (cmf)
```


### Child, mother, and father trio with parent-of-origin effect

The `cmf` model can be extended to capture a parent-of-origin effect from the mother or the father, `cmf_mt` and `cmf_ft`, respectively.

```
y = �m m + �c c + �f f + �mt hmt + e                                                                  (cmf_mt)
```
```
y = �m m + �c c + �f f + �ft hft + e                                                                  (cmf_mt)
```


### Sub-models

From these models, the following can be obtained by removing degrees of freedom.


```
y = �m m + �c c + e                                                                                   (cm)
```

```
y = �m m + �c c + �mt hmt + e                                                                         (cm_mt)
```

```
y = �m m + �c c + �ft hft + e                                                                         (cm_ft)
```

```
y = �c c + �f f + e                                                                                   (cf)
```

```
y = �c c + �f f + �mt hmt + e                                                                         (cf_mt)
```

```
y = �c c + �f f + �ft hft + e                                                                         (cf_ft)
```

```
y = �m m + �f f + e                                                                                   (mf)
```

```
y = �c c + e                                                                                          (c)
```

```
y = �c c + �mt hmt + e                                                                                (c_mt)
```

```
y = �c c + �ft hft + e                                                                                (c_ft)
```

```
y = �m m + e                                                                                          (m)
```

```
y = �m m + �mt hmt + e                                                                                (m_mt)
```

```
y = �f f + e                                                                                          (f)
```

```
y = �f f + �ft hft + e                                                                                (f_mt)
```

