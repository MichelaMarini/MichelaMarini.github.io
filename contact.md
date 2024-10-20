---
layout: page
---

<!-- No Title -->

<style>
  body {
    background-color: #E1EFC1; 
  }

  a {
    color: #00008B; /* Dark blue for links */
    background-color: transparent; /* Ensure link background is transparent */
    text-decoration: none; /* Remove underline, optional */
  }

  a:visited {
    color: #5A5A8B; /* Slightly lighter dark blue for visited links */
  }

  .contact-container {
    display: flex;
    justify-content: flex-start; /* Align content to the left */
    align-items: center;
    height: 50vh;
  }

  .contact-info {
    max-width: 400px;
    text-align: left;
  }

 .contact-image-container {
    float: left; /* Float the image container to the left on larger screens */
    margin-left: 80px; /* Slightly reduce the left margin */
    border-radius: 50%; /* Keep the image circular */
    width: 240px; /* Reduce the width to make the image smaller */
    height: 240px; /* Reduce the height accordingly */
    overflow: hidden; /* Hide the overflow to maintain the circular shape */
}

  .contact-image {
    width: 100%; /* Ensure the image fits within the container */
    height: auto; /* Allow the image to adjust its height while maintaining aspect ratio */
  }

/* Responsive adjustments */
@media (max-width: 768px) {
    .contact-image-container {
        float: none; /* Remove float on smaller screens */
        margin-left: 0; /* Remove the left margin */
        width: 110px; /* Further reduce the width for smaller screens */
        height: 110px; /* Further reduce the height for smaller screens */
        margin: 0 auto; /* Center the image */
    }

    .contact-image {
        width: 100%; /* Ensure the image fits the container */
        height: 100%; /* Ensure the image fills the container */
        object-fit: cover; /* Keep the circular aspect ratio */
        display: block;
    }
}

<style>
  .contact-container .contact-info h1 {
    font-family: 'Arial', sans-serif; /* Change this to any font you prefer */
    font-size: 24px; /* Adjust the size */
    font-weight: bold; /* Adjust the weight */
    color: #333; /* Adjust the color */
    border-radius: 50%; /* Ensure circular shape */
  }
</style>

<div class="contact-container">
  <div class="contact-info">
    <h1>Get in Touch!</h1>

    <ul>
      <li><b>Email:</b> <a href="mailto:mmarini2@uh.edu">mmarini2(at)uh(dot)edu</a></li>
      <li><b>LinkedIn:</b> <a href="https://www.https://www.linkedin.com/in/michela-marini-bb813517a/">Michela Marini</a></li>
      <li><b>GitHub:</b> <a href="https://github.com/MichelaMarini" target="_blank">MichelaMarini</a></li>
      <li><b>Office:</b> <a href="https://www.uh.edu/maps/#/find/PGH">PGH</a> 612</li>
      <li><b>Address:</b> 3551 Cullen Blvd., Houston, TX 77004</li>
    </ul>
  </div>

  <div class="contact-image-container">
    <img src="/images/me.jpeg" alt="Contact Image" class="contact-image">
  </div>
</div>
