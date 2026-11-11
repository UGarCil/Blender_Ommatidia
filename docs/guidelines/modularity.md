When writing programs, it's important to ensure proper separation fo concerns and a hierarchical structure that encapsulates testable behavior. In code_practices.md we establish the practices that ensure scalable code by using the How To Design Recipes principle, requiring a .md that shows a clear understanding on the information domain in the construction of atomic algorithms, and how the scale to create modules in what I called compound algorithms.

This documents outlines the best practices to create the software infrastructure: how to organize the folders, how to integrate testing based on the data definitions.

In scientific pipelines we encounter an unusual limitation while writing good programs: we don't completely understand the full breath on the information domain! The results of the analysis will be used in a positive loop to reshape the base code. So it's important to build for resilience, and I like the analogy of lego blocks: independent, clearly defined units of information that can be stacked and assembles, properly evaluated independently and replaced easily. 


### The role of serialization and data persistence
When dealing with large data corpora we deal with the need of optimization, parallelization and the sort, drawing upon the writing of persistent files to ensure data is not lost. All of these considerations open the gate to additional challenges:

- The HtdP templates aid in the identification breaking points for parallelization and optimization: The clear use of data definitions and modularization ensure that the atomic unit of parallelization can be promptly found and targeted at the right level of modularity (i.e. which modules upstream in the chain of compound algorithms are affected). If the program was built right, then the division where to cut should be obvious, often at the module.py that is affected. This will require assessment from agents with human assistance involved to ensure a compound understanding. This is the best way to achieve a goal.
- building for scale means to ensure the serialization is independent, and work is resumable if a node or thread is lost. For this, we'll need to have additional modules dedicated to the validation and health of data stored. 


### Staying organized
Because of the sequential nature of the data transformations, it is important to introduce a sense of serialization on the pipeline. While sometimes I like using enumerations like 00_import.py, 01_read.py, 02_plot.py, this is not a desirable behavior when trying to build a package. Instead, we will delegate this sequence to **the order in which we write statements in the main.py orchestrator file**. Such main.py file will live outside of the package src, at the root, and will be used to call the pipeline. The main.py will be fully explained in the documentation as well, which we will build inside docs/user. Such documentation will have a user-focused guide rather than details on the modularization of the program (although of course sometimes it may be necessary to reference).

The organization of our project follows the standard src/ folder under which we write python packages, accompanied by a pyproject.toml. As follows:

```
root/
    src/
        package_name/
            compound_algorithm_module_name/
                module_name.py
                module_name.md
                atomic_algorithm_scripts/
                    script_name1.py
                    script_name2.py
    main.py (the only script exposed at root and orchestrator module)    
    /test
```

main.py is designed to be simple and easy to read. This will be a guiding system for the end user to know how to use our software. 

Keep in mind, the absolute path of the folders is subject to change, as I may want to move the folder root around. However, the structure inside of root will not change, so we should keep relative paths in mind.

When traversing or walking our paths, we can rely on the principles dictated by our templates as well.

We also need to clearly specifying where the input and output directories go. In the use of main.py, it will be desirable to run the results from analyses, simulations, and essentially the running of the program in a different folder. Canonically I refer to such folder as 00_results (the digits place it at the top of the folder list). 


### Documentation

Each module will have a .md referring to it, and in the folder where such module .py is. This means that connections between our module documentation will be stated via relative paths (and opening in a new document rather than replacing the text where the user clicks, and equivalent to _blank in HTML pages), allowing the user to stay in a folder where the code lives. In addition, such modularity will be presented to the end user as details for the package under docs/user/modules.md. Such document will be an introductory entrypoint to the modularity of the program and will include a diagram and an easy to understand disclosure on such modularity. 




### Testing

Each module has to be tested. We will implement thorough testing based on the types and edge cases that come from the data definitions we use. We will begin at the bottom and at an atomic level. First the atomic data definitions, where we can apply mutational and edge-case tests to determine what breaks the program. We will be using our recipes to implement clear validators in the codebase when bugs are found.

### Keeping a logbook

This goes under agent/ and it's for agents to keep a log (human readable) on our findings, highlights, and decisions made along the way. Such logbook will expand as the proyect moves forward.

### Writing your documents inside docs/agents/

The files docs/agent/data_definitions.md and modules.md are intially left empty and were created by the user. If not already, fill them up with the information you need to record a thorough analysis on the actual solution to the problem (i.e. actualy write the data definitions, modular structure, and anything else). If the files are written already you can:

- Use them as a guideline. 
- Propose your human support for feedback if you need to change the atomic units (should happen rarely)
- Don't replace the guides without consulting

Why? In the atomic units and modularity we have the axioms we'll follow to solve a problem. Changing them involves changing everything upstream, and should be a last resort once our understanding of the problem reaches convergence (i.e. you and I agreed on it).

### Writing a manuscript in latex

In the root you'll also find a folder called manuscript/

That manuscript provides an ecosystem with everything needed to write a scientific-writing styled paper with the results and findings of our program. This is important as a lot of the information found in the booklog will end up here, presented for a scientific audience. Inside manuscript we shoudl feel free to create a .venv if needed called .venv_man dedicated exclusively to work on the manuscript. Inside ./manuscript/scripts we will write scripts properly named and allocated inside of a folder structure as follows:

```
./manuscript
    /scripts
        plot_metric_name1/
            plot_metric_name1.py
            results/
        plot_metric_name2/
            plot_metric_name2.py
            results/
        
```