// This script fixes the Run Pipelines tab content issue

document.addEventListener('DOMContentLoaded', function() {
  // Function to update Run Pipelines tab content
  function updateRunPipelinesContent() {
    const configContent = document.getElementById('config-content');
    if (!configContent) return;
    
    // Check if Demo environment is selected
    const selectedEnv = document.querySelector('#environment-options .option.selected');
    if (!selectedEnv) return;
    
    const isDemoSelected = selectedEnv.getAttribute('data-value') === 'demo';
    
    // If not Demo, update content
    if (!isDemoSelected) {
      // Save the heading
      const heading = configContent.querySelector('h2');
      const headingText = heading ? heading.textContent : 'Run Pipelines';
      
      // Update content immediately
      configContent.innerHTML = `
        <h2>${headingText}</h2>
        
        <div class="tracer-instructions" style="margin-top: 2rem; border: 1px solid var(--border); padding: 1.5rem; border-radius: 8px; background-color: rgba(0,0,0,0.2);">
          <h4>Step 3: Run pipeline</h4>
          <p>We have pre-installed some pipelines the github repo for you to run.</p>
          <p>We would recommend to start with a simple rnaseq pipeline in Nextflow:</p>
          <div class="code-container">
            <pre><code>nextflow run nf-core/rnaseq -c custom.config -profile docker,test --outdir</code></pre>
            <button class="copy-btn">
              <i class="fas fa-copy"></i>
            </button>
          </div>
          <p>Other pipelines, written in Bash, Nextflow, WDL, and CWL can be found under the <code>pipelines</code> files for you to test out.</p>
          <p>Play around with the other pipelines, have fun! 🧬</p>
        </div>
      `;
      
      // Add event listeners to copy buttons
      configContent.querySelectorAll('.copy-btn').forEach(btn => {
        btn.addEventListener('click', function() {
          const code = this.previousElementSibling.textContent;
          navigator.clipboard.writeText(code);
          
          // Visual feedback
          const originalIcon = this.innerHTML;
          this.innerHTML = '<i class="fas fa-check"></i>';
          setTimeout(() => {
            this.innerHTML = originalIcon;
          }, 2000);
        });
      });
    }
  }
  
  // Completely override the tab click behavior for the Run Pipelines tab
  const configTab = document.querySelector('.tab[data-tab="config"]');
  if (configTab) {
    configTab.addEventListener('click', function(event) {
      // Get all tabs and remove active class
      document.querySelectorAll('.tab').forEach(tab => {
        tab.classList.remove('active');
      });
      
      // Add active class to this tab
      this.classList.add('active');
      
      // Hide all tab content
      document.querySelectorAll('.tab-content').forEach(content => {
        content.classList.remove('active');
      });
      
      // Get the config content element
      const configContent = document.getElementById('config-content');
      
      // Show the config content
      configContent.classList.add('active');
      
      // Check if Demo environment is selected
      const selectedEnv = document.querySelector('#environment-options .option.selected');
      if (!selectedEnv) return;
      
      const isDemoSelected = selectedEnv.getAttribute('data-value') === 'demo';
      
      // If not Demo, update content immediately
      if (!isDemoSelected) {
        // Save the heading or create a new one
        const headingText = 'Run Pipelines';
        
        // Update content immediately
        configContent.innerHTML = `
          <h2>${headingText}</h2>
          
          <div class="tracer-instructions" style="margin-top: 2rem; border: 1px solid var(--border); padding: 1.5rem; border-radius: 8px; background-color: rgba(0,0,0,0.2);">
            <h4>Step 3: Run pipeline</h4>
            <p>We have pre-installed some pipelines the github repo for you to run.</p>
            <p>We would recommend to start with a simple rnaseq pipeline in Nextflow:</p>
            <div class="code-container">
              <pre><code>nextflow run nf-core/rnaseq -c custom.config -profile docker,test --outdir</code></pre>
              <button class="copy-btn">
                <i class="fas fa-copy"></i>
              </button>
            </div>
            <p>Other pipelines, written in Bash, Nextflow, WDL, and CWL can be found under the <code>pipelines</code> files for you to test out.</p>
            <p>Play around with the other pipelines, have fun! 🧬</p>
          </div>
        `;
        
        // Add event listeners to copy buttons
        configContent.querySelectorAll('.copy-btn').forEach(btn => {
          btn.addEventListener('click', function() {
            const code = this.previousElementSibling.textContent;
            navigator.clipboard.writeText(code);
            
            // Visual feedback
            const originalIcon = this.innerHTML;
            this.innerHTML = '<i class="fas fa-check"></i>';
            setTimeout(() => {
              this.innerHTML = originalIcon;
            }, 2000);
          });
        });
      } else {
        // For Demo environment, let the original code handle it
        // This will show the original Demo content
      }
      
      // Prevent the default click behavior
      event.preventDefault();
      event.stopPropagation();
    }, true); // Use capture phase to intercept the event before other handlers
  }
  
  // Also run when environment options are clicked
  document.querySelectorAll('#environment-options .option').forEach(option => {
    option.addEventListener('click', function() {
      // Check if Run Pipelines tab is active
      const isConfigTabActive = document.querySelector('.tab[data-tab="config"]').classList.contains('active');
      if (isConfigTabActive) {
        // Update content immediately
        updateRunPipelinesContent();
      }
    });
  });
});

// Add this function to fix the Install Tracer tab content
document.addEventListener('DOMContentLoaded', function() {
  // Function to fix the Install Tracer tab content
  function fixInstallTracerContent() {
    const installContent = document.getElementById('install-content');
    if (!installContent) return;
    
    // Check if Demo environment is selected
    const selectedEnv = document.querySelector('#environment-options .option.selected');
    if (!selectedEnv) return;
    
    const isDemoSelected = selectedEnv.getAttribute('data-value') === 'demo';
    
    // Check if "Run Your Own" is selected in Step 4
    const runYourOwnSelected = document.querySelector('#pipeline-options .option[data-value="own"].selected') !== null;
    
    // If not Demo, update content to remove Step 3
    if (!isDemoSelected) {
      // Find all h4 elements
      const h4Elements = installContent.querySelectorAll('h4');
      
      // Find the Step 3 heading
      let step3Heading = null;
      for (const h4 of h4Elements) {
        if (h4.textContent.includes('Step 3: Run pipeline')) {
          step3Heading = h4;
          break;
        }
      }
      
      // If Step 3 heading is found, remove it and all elements until the next heading or end
      if (step3Heading) {
        // Find the parent container of Step 3
        let parentContainer = step3Heading.parentElement;
        
        // Remove Step 3 section from the parent container
        parentContainer.innerHTML = parentContainer.innerHTML.replace(
          /<h4[^>]*>Step 3: Run pipeline<\/h4>[\s\S]*?(?=<h4|<\/div>)/g, 
          ''
        );
      }
      
      // Also remove the text about other pipelines and playing around
      const paragraphs = installContent.querySelectorAll('p');
      for (const p of paragraphs) {
        if (p.textContent.includes('Other pipelines, written in Bash, Nextflow, WDL, and CWL can be found under the') ||
            p.textContent.includes('Play around with the other pipelines, have fun!')) {
          p.remove();
        }
        
        // Replace [exclamation mark] with ❗
        if (p.textContent.includes('[exclamation mark]')) {
          p.innerHTML = p.innerHTML.replace('[exclamation mark]', '❗');
        }
      }
      
      // If "Run Your Own" is selected, completely replace the content
      if (runYourOwnSelected) {
        // Find the tracer instructions container
        const tracerInstructions = installContent.querySelector('.tracer-instructions');
        
        if (tracerInstructions) {
          // Save the heading if it exists
          const heading = installContent.querySelector('h2');
          const headingText = heading ? heading.textContent : 'Install Tracer';
          
          // Replace the entire content with our custom content
          tracerInstructions.innerHTML = `
            <p>After installing the tracer daemon, you can run any pipeline you want</p>
            <p>One thing to add before the start of a run is the below line:</p>
            <div class="code-container">
              <pre><code>tracer init</code></pre>
              <button class="copy-btn">
                <i class="fas fa-copy"></i>
              </button>
            </div>
            <p>It will ask you to setup the parameters of this pipeline. This allows for customisation and easy information retrieval later:</p>
            <ul style="list-style-type: none; padding-left: 1rem;">
              <li>--pipeline-name: what name do you give this pipeline? (one pipeline consists of multiple runs)</li>
              <li>--environment: what is the environment you are running in?</li>
              <li>--pipeline-type: what type of pipeline are you running e.g., rnaseq, scrnaseq, xxxx.<br>
                  protip: you could also name this after a cost centre or team so that you can easily search the other pipelines</li>
              <li>--user-operator: what is your name/who is running this pipeline?</li>
            </ul>
            <p>Tip: add this line to the start of your pipeline file. This way, you never forget to start monitoring your pipeline</p>
          `;
          
          // Add event listeners to copy buttons
          tracerInstructions.querySelectorAll('.copy-btn').forEach(btn => {
            btn.addEventListener('click', function() {
              const code = this.previousElementSibling.textContent;
              navigator.clipboard.writeText(code);
              
              // Visual feedback
              const originalIcon = this.innerHTML;
              this.innerHTML = '<i class="fas fa-check"></i>';
              setTimeout(() => {
                this.innerHTML = originalIcon;
              }, 2000);
            });
          });
        }
      }
    }
  }
  
  // Run the function when environment options are clicked
  document.querySelectorAll('#environment-options .option').forEach(option => {
    option.addEventListener('click', function() {
      // Fix the Install Tracer tab content
      setTimeout(fixInstallTracerContent, 100);
    });
  });
  
  // Run the function on page load
  setTimeout(fixInstallTracerContent, 100);
  
  // Add event listener for the Install Tracer tab
  const installTab = document.querySelector('.tab[data-tab="install"]');
  if (installTab) {
    installTab.addEventListener('click', function() {
      // Fix the Install Tracer tab content when the tab is clicked
      setTimeout(fixInstallTracerContent, 100);
    });
  }
  
  // Add event listeners to pipeline options
  document.querySelectorAll('#pipeline-options .option').forEach(option => {
    option.addEventListener('click', function() {
      // Fix the Install Tracer tab content when pipeline options are clicked
      setTimeout(fixInstallTracerContent, 100);
    });
  });
  
  // Override the ensureInstallTracerContent function if it exists
  if (typeof window.ensureInstallTracerContent === 'function') {
    const originalEnsureInstallTracerContent = window.ensureInstallTracerContent;
    
    window.ensureInstallTracerContent = function() {
      // Call the original function
      originalEnsureInstallTracerContent();
      
      // Then apply our fix
      setTimeout(fixInstallTracerContent, 10);
    };
  }
});
