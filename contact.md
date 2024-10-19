---
layout: page
---

<!-- No Title -->

<style>
  body {
    background-color: #add8e6; /* Light blue background */
  }

  a {
    color: #00008B; /* Dark blue for links */
    background-color: transparent; 
    text-decoration: none; 
  }

  a:visited {
    color: #5A5A8B; 
  }
  </style>

  /* Optional: if links are inside any specific containers that might have background issues */
  .container, .content, .wrapper {
    background-color: transparent; /* Ensure parent containers have transparent backgrounds */
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
    float: left; /* Float the image container to the left */
    margin-left: 100px; /* Add space between the contact info and the image */
    border-radius: 50%; /* Make the image circular */
    width: 300px; /* Adjust the width to your preference */
    height: 300px; /* Adjust the height to your preference */
    overflow: hidden; /* Hide the overflow to maintain the circular shape */
  }

  .contact-image {
    width: 100%; /* Ensure the image fits within the container */
    height: auto; /* Allow the image to adjust its height while maintaining aspect ratio */
  }

<style>
  .contact-container .contact-info h1 {
    font-family: 'Arial', sans-serif; /* Change this to any font you prefer */
    font-size: 24px; /* Adjust the size */
    font-weight: bold; /* Adjust the weight */
    color: #333; /* Adjust the color */
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
