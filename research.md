---
title: Research
---

My research focuses on developing computational methods to process and analyze biological images. I combine applied mathematics, statistical analysis, and machine learning to study complex biological data, with a particular interest in microscopy imaging and spatial omics. Through my work, I aim to develop tools that help researchers better understand cellular features, biological phenotypes, and tissue organization.

<!-- Edit the three themes below to describe your work in your own words. -->
<div class="themes">
  <article class="theme">
    <h3>Image analysis and machine learning</h3>
    <p>Supervised and unsupervised learning methods to quantify cell morphology and structural changes from microscopy images.</p>
  </article>
  <article class="theme">
    <h3>Single-cell phenotyping</h3>
    <p>High-throughput, image-based profiling of Cell Painting data, including open-source software for single-cell analysis.</p>
  </article>
  <article class="theme">
    <h3>Spatial biology</h3>
    <p>Computational methods for spatial omics data to study how cells are organized within tissues.</p>
  </article>
</div>

[Interactive analysis of the Cell Painting dataset](/cpmetab/)

## Recent publications

<ol class="pubs">
{% for pub in site.data.publications limit: 3 %}{% include publication.html pub=pub %}{% endfor %}
</ol>

<p><a class="more" href="{{ '/publications/' | relative_url }}">See all publications</a></p>
