---
title: Contact
---

<section class="contact">
  <img src="{{ '/images/me.jpeg' | relative_url }}" alt="Portrait of {{ site.author.name }}" class="portrait portrait-round">
  <div>
    <p class="lead">Feel free to get in touch about research, collaborations or software.</p>
    <dl class="facts">
      <dt>Email</dt><dd><a href="mailto:{{ site.author.email }}">{{ site.author.email }}</a></dd>
      <dt>Group</dt><dd><a href="{{ site.author.group_url }}">{{ site.author.group }}</a></dd>
      <dt>Institute</dt><dd>{{ site.author.institute }}</dd>
      <dt>City</dt><dd>{{ site.author.city }}</dd>
    </dl>
    {% include social-links.html %}
  </div>
</section>
