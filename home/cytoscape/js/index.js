var x = [1,2,3,4,5,6,7,8,9];
var y = [0,0,0,0,0,0,0,0,0];

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

    let $color = $('#color');
    let getStylesheet = name => {
      let convert = res => toJson(res);

      return fetch(`stylesheets/base.json`).then( convert );
    };
    let applyStylesheet = stylesheet => {

      cy.style().fromJson( stylesheet ).update();

      if ($color.value == "stage"){
        var y = "mapData("+$color.value+",6.5,8.5,blue,red)"
      }else{
        var y = "mapData("+$color.value+",0,1,blue,red)"
      }
      cy.style().selector("node").style({'background-color':y}).update()
    };
    let applyStylesheetFromSelect = () => Promise.resolve( $color.value ).then( getStylesheet ).then( applyStylesheet );

    let $weights = $('#weights');
    let getWeights = name => {
      let convert = res => toJson(res);

      return fetch(`stylesheets/base.json`).then( convert );
    };
    let applyWeights = weight => {

      cy.style().fromJson( weight ).update();

      switch ($weights.value){
        case "weight": var y = "mapData(weight,0,1,blue,red)"; break;
        case "weight2": var y = 3; break;
      }
      cy.style().selector("edge").style({'width':y}).update()
    };
    let applyWeightsFromSelect = () => Promise.resolve( $weights.value ).then( getWeights ).then( applyWeights );

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
    let applyLayoutFromSelect = () => Promise.resolve( "hierarchical" ).then( getLayout ).then( applyLayout );

    cy = window.cy = cytoscape({
      container: $('#cy'),

      layout: { // replace with your own layout parameters
          name: 'breadthfirst',
          directed: true,
          padding: 30
      }
    });

    tryPromise( applyDatasetFromSelect ).then( applyStylesheetFromSelect ).then( applyWeightsFromSelect ).then( applyLayoutFromSelect );

    $dataset.addEventListener('change', function(){
      tryPromise( applyDatasetFromSelect ).then( applyLayoutFromSelect );
    });

    $color.addEventListener('change', applyStylesheetFromSelect);

    $weights.addEventListener('change', applyWeightsFromSelect);

    // $('#redo-layout').addEventListener('click', applyLayoutFromSelect);

    let $a = $('#hola');
    cy.bind('click', 'node', function(node) {
      Promise.resolve($a)
      // console.log(node.cyTarget.predecessors().edges());
      switch($a.value){
        case "false":

          break;
        case "backward":
          node.cyTarget.incomers().edges().animate({
            style: {
              lineColor: "red"
            }
          });
          break;
        case "forward":
          node.cyTarget.outgoers().edges().animate({
            style: {
              lineColor: "red"
            }
          });
          break;
        case "expression":
          Promise.resolve($color)
          x 
          y[] = node.cyTarget.data($color.value)
          console.log(node.cyTarget.data($color.value))
          Plotly.newPlot('myDiv', plotData);
          break;
        }
      let name = node.cyTarget.data("id");
      document.getElementById('annotation').innerHTML = "<b>Node: </b>"+name;
      console.log(name)
    });
    cy.bind('cxttap', 'node', function(node) {
      // console.log(node.cyTarget.predecessors().edges());
      node.cyTarget.incomers().edges().animate({
        style: {
          lineColor: "lightgrey"
        }
      });
      console.log(node.cyTarget.predecessors().edges());
      node.cyTarget.outgoers().edges().animate({
        style: {
          lineColor: "lightgrey"
        }
      });
    });

    // cy.on("click", "node", (evt) => {evt.cyTarget.ancestors().animate({
    //   style: {lineColor: "red"} }) })
    // cy.off("tap", "node", (evt) => {evt.cyTarget.successors().animate({
    //     style: {lineColor: "grey"} }) })
  
    var plotData = [
      {
        x: x,
        y: y,
        error_y: {
          type: 'data',
          array: [1, 2, 3],
          visible: true
        },
        type: 'scatter'
      }
    ];
    console.log(plotData[0]["x"][1])
    Plotly.newPlot('myDiv', plotData);

  });
  
})();

// tooltips with jQuery
$(document).ready(() => $('.tooltip').tooltipster());