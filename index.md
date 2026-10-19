---
title: Home
---

<section class="intro">
  <div class="intro-text" markdown="1">

I am a Postdoctoral Researcher in the [Spatial Biotechnology group](https://ibecbarcelona.eu/research-groups/spatial-biotechnology/) at the Institute for Bioengineering of Catalonia (IBEC), working at the intersection of computational imaging, spatial biology, and applied mathematics. I hold a Ph.D. in Applied Mathematics from the University of Houston, USA.

My research focuses on developing computational methods and software to process and analyze biological images and spatial omics data. I enjoy using image processing, machine learning, and quantitative methods to solve challenging problems in biology and better understand cellular features, phenotypes, and tissue organization.

Previously, I collaborated with the [Integrated Microscopy Core Lab](https://www.bcm.edu/research/atc-core-labs/integrated-microscopy-core) at Baylor College of Medicine on the quantitative analysis of Cell Painting images. I also served as vice president of the University of Houston chapter of the Association for Women in Mathematics (AWM).

  </div>
  <aside class="intro-card">
    <img src="{{ '/images/me.jpeg' | relative_url }}" alt="Portrait of {{ site.author.name }}" class="portrait">
    <dl class="facts">
      <dt>Position</dt><dd>Postdoc</dd>
      <dt>Group</dt><dd><a href="{{ site.author.group_url }}">{{ site.author.group }}</a></dd>
      <dt>Institute</dt><dd>IBEC, {{ site.author.city }}</dd>
      <dt>Email</dt><dd><a href="mailto:{{ site.author.email }}">{{ site.author.email }}</a></dd>
    </dl>
    {% include social-links.html %}
  </aside>
</section>

## News

<ul class="news">
{% assign news = site.data.news | sort: "date" | reverse %}
{% for item in news %}
  <li><time datetime="{{ item.date | date: '%Y-%m-%d' }}">{{ item.date | date: "%b %Y" }}</time><span>{{ item.text }}</span></li>
{% endfor %}
</ul>

## Education

<ul class="timeline">
  <li><span class="degree">Ph.D. in Mathematics</span><span class="place">University of Houston, USA</span></li>
  <li><span class="degree">M.S. in Bioelectronics Engineering</span><span class="place">University of Bologna, Italy</span></li>
  <li><span class="degree">B.S. in Biomedical Engineering</span><span class="place">University of Bologna, Italy</span></li>
</ul>
