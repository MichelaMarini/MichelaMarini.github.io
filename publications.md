---
title: Publications
---

{% assign by_year = site.data.publications | group_by: "year" | sort: "name" | reverse %}
{% for group in by_year %}
<h2 class="year">{{ group.name }}</h2>
<ol class="pubs">
{% for pub in group.items %}{% include publication.html pub=pub %}{% endfor %}
</ol>
{% endfor %}
