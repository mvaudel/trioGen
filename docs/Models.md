## Trio Models

TrioGen proposes a selection of {child, mother, father} trio models for trio and haplotype analyses. Please note that it is currently not possible to implement custom models. If this is something important for your research, please [open an issue](https://github.com/mvaudel/trioGen/issues).


### Haplotypes

In the haplotype `h` model, the regression is conducted against the number of transmitted and non-transmitted tested alleles, h, as done by [Chen _et al._](https://doi.org/10.1101/737106).

```
y = ßmnt hmnt + ßmt hmt + ßft hft + ßfnt hfnt + e                                                     (h)
```

Where y represents the phenotypes, _hmnt_ the number of maternal tested alleles non-transmitted to the child, _hmt_ the number of maternal tested alleles transmitted to the child, _hft_ the number of paternal tested alleles transmitted to the child, and _hfnt_ the number of paternal tested alleles non-transmitted to the child.


### Child, mother, and father trio

In the {child, mother, father} `cmf` model, the regression is conducted against the number of tested alleles in the child, mother, and father conditioned against each other.

```
y = ßm m + ßc c + ßf f + e                                                                            (cmf)
```


### Child, mother, and father trio with parent-of-origin effect

The `cmf` model can be extended to capture a parent-of-origin effect from the mother or the father, `cmf_mt` and `cmf_ft`, respectively.

```
y = ßm m + ßc c + ßf f + ßmt hmt + e                                                                  (cmf_mt)
```
```
y = ßm m + ßc c + ßf f + ßft hft + e                                                                  (cmf_mt)
```


### Sub-models

From these models, the following can be obtained by removing degrees of freedom.


```
y = ßm m + ßc c + e                                                                                   (cm)
```

```
y = ßm m + ßc c + ßmt hmt + e                                                                         (cm_mt)
```

```
y = ßm m + ßc c + ßft hft + e                                                                         (cm_ft)
```

```
y = ßc c + ßf f + e                                                                                   (cf)
```

```
y = ßc c + ßf f + ßmt hmt + e                                                                         (cf_mt)
```

```
y = ßc c + ßf f + ßft hft + e                                                                         (cf_ft)
```

```
y = ßm m + ßf f + e                                                                                   (mf)
```

```
y = ßc c + e                                                                                          (c)
```

```
y = ßc c + ßmt hmt + e                                                                                (c_mt)
```

```
y = ßc c + ßft hft + e                                                                                (c_ft)
```

```
y = ßm m + e                                                                                          (m)
```

```
y = ßm m + ßmt hmt + e                                                                                (m_mt)
```

```
y = ßf f + e                                                                                          (f)
```

```
y = ßf f + ßft hft + e                                                                                (f_mt)
```

