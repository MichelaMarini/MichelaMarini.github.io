---
title: Contact
---

<section class="contact">
  <img src="{{ '/images/me.jpeg' | relative_url }}" alt="Portrait of {{ site.author.name }}" class="portrait portrait-round">
  <div>
    <dl class="facts">
      <dt>Email</dt><dd><a href="mailto:{{ site.author.email }}">{{ site.author.email }}</a></dd>
      <dt>City</dt><dd>{{ site.author.city }}</dd>
    </dl>
    {% include social-links.html %}
  </div>
</section>
