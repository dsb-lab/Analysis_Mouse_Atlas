(function(){
  document.addEventListener('DOMContentLoaded', function(){
    let $$ = selector => Array.from( document.querySelectorAll( selector ) );
    let $ = selector => document.querySelector( selector );

    let tryPromise = fn => Promise.resolve().then( fn );

    let toJson = obj => obj.json();
    let toText = obj => obj.text();

    let cy;

    let $dataset = $('#data');
    let getDataset = name => fetch(`datasets/${name}`).then( toJson );
    let applyDataset = dataset => {
      // so new eles are offscreen
      cy.zoom(0.001);
      cy.pan({ x: -9999999, y: -9999999 });

      // replace eles
      cy.elements().remove();
      cy.add( dataset );
    }
    let applyDatasetFromSelect = () => Promise.resolve( $dataset.value ).then( getDataset ).then( applyDataset );

    let $stylesheet = $('#color');
    let getStylesheet = name => {
      let convert = res => toJson(res);

      return fetch(`stylesheets/stage.json`).then( convert );
    };
    let applyStylesheet = stylesheet => {
      if( typeof stylesheet === typeof '' ){
        cy.style().fromString( stylesheet ).update();
      } else {
        cy.style().fromJson( stylesheet ).update();
      }
      if ($stylesheet.value == "stage"){
        var y = "mapData("+$stylesheet.value+",6.5,8.5,blue,red)"
      }else{
        var y = "mapData("+$stylesheet.value+",0,1,blue,red)"
      }
      cy.style().selector("node").style({'background-color':y}).update()
    };
    let applyStylesheetFromSelect = () => Promise.resolve( $stylesheet.value ).then( getStylesheet ).then( applyStylesheet );

    let $layout = $('#layout');
    let layouts = {
      hierarchical: { // replace with your own layout parameters
        name: 'breadthfirst',
        directed: true,
        padding: 30
      }
    };
    let prevLayout;
    let getLayout = name => Promise.resolve( layouts[ name ] );
    let applyLayout = layout => {
      if( prevLayout ){
        prevLayout.stop();
      }

      let l = prevLayout = cy.makeLayout( layout );

      return l.run().promiseOn('layoutstop');
    }
    let applyLayoutFromSelect = () => Promise.resolve( $layout.value ).then( getLayout ).then( applyLayout );

    let $algorithm = $('#algorithm');
    let getAlgorithm = (name) => {
      switch (name) {
        case 'bfs': return Promise.resolve(cy.elements().bfs.bind(cy.elements()));
        case 'dfs': return Promise.resolve(cy.elements().dfs.bind(cy.elements()));
        case 'astar': return Promise.resolve(cy.elements().aStar.bind(cy.elements()));
        case 'none': return Promise.resolve(undefined);
        case 'custom': return Promise.resolve(undefined); // replace with algorithm of choice
        default: return Promise.resolve(undefined);
      }
    };
    let runAlgorithm = (algorithm) => {
      if (algorithm === undefined) {
        return Promise.resolve(undefined);
      } else {
        let options = {
          root: '#' + cy.nodes()[0].id(),
          // astar requires target; goal property is ignored for other algorithms
          goal: '#' + cy.nodes()[Math.round(Math.random() * (cy.nodes().size() - 1))].id()
        };
        return Promise.resolve(algorithm(options));
      }
    }
    let currentAlgorithm;
    let animateAlgorithm = (algResults) => {
      // clear old algorithm results
      cy.$().removeClass('highlighted start end');
      currentAlgorithm = algResults;
      if (algResults === undefined || algResults.path === undefined) {
        return Promise.resolve();
      }
      else {
        let i = 0;
        // for astar, highlight first and final before showing path
        if (algResults.distance) {
          // Among DFS, BFS, A*, only A* will have the distance property defined
          algResults.path.first().addClass('highlighted start');
          algResults.path.last().addClass('highlighted end');
          // i is not advanced to 1, so start node is effectively highlighted twice.
          // this is intentional; creates a short pause between highlighting ends and highlighting the path
        }
        return new Promise(resolve => {
          let highlightNext = () => {
            if (currentAlgorithm === algResults && i < algResults.path.length) {
              algResults.path[i].addClass('highlighted');
              i++;
              setTimeout(highlightNext, 500);
            } else {
              // resolve when finished or when a new algorithm has started visualization
              resolve();
            }
          }
          highlightNext();
        });
      }
    };
    let applyAlgorithmFromSelect = () => Promise.resolve( $algorithm.value ).then( getAlgorithm ).then( runAlgorithm ).then( animateAlgorithm );

    cy = window.cy = cytoscape({
      container: $('#cy')
    });

    tryPromise( applyDatasetFromSelect ).then( applyStylesheetFromSelect ).then( applyLayoutFromSelect );

    $dataset.addEventListener('change', function(){
      tryPromise( applyDatasetFromSelect ).then( applyLayoutFromSelect ).then ( applyAlgorithmFromSelect );
    });

    $stylesheet.addEventListener('change', applyStylesheetFromSelect);

    $layout.addEventListener('change', applyLayoutFromSelect);

    $('#redo-layout').addEventListener('click', applyLayoutFromSelect);

    $algorithm.addEventListener('change', applyAlgorithmFromSelect);

    $('#redo-algorithm').addEventListener('click', applyAlgorithmFromSelect);
  });
})();

// tooltips with jQuery
$(document).ready(() => $('.tooltip').tooltipster());