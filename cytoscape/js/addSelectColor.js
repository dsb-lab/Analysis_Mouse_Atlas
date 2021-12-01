//Codigo a Ejecutar al Cargar la Pagina
function myOnLoad() {
    loadColorCodes()
   }
   
// funcion para Cargar Provincias al campo <select>
function loadColorCodes() {

    var array = [
    "stage",
    'Adh1a2',
    'Bmp4',
    'Bmp8',
    'Cdh1',
    'Cdh2',
    'Cdh11',
    'Cer1',
    'Dlx5',
    'Dmrt2',
    'Ebf1',
    'Egr2',
    'En1',
    'Eomes',
    'Epcam',
    'Eya1',
    'Eya2',
    'Evx1',
    'Fgf5',
    'Foxa2',
    'Foxc1',
    'Foxd1',
    'Foxf1',
    'Foxf2',
    'Foxp2',
    'Hand1',
    'Hand2',
    'Hesx1',
    'Hoxa1',
    'Hoxa2',
    'Hoxa3',
    'Hoxa4',
    'Hoxa10',
    'Irx3',
    'Isl1',
    'Kdr',
    'Lefty1',
    'Lef1',
    'Lfng',
    'Lhx1',
    'Meox1',
    'Mesp1',
    'Mesp2',
    'Msc',
    'Msgn',
    'Myf5',
    'Nanog',
    'Nkx2-5',
    'Noto',
    'Olig2',
    'Osr1',
    'Otx2',
    'Pax1',
    'Pax3',
    'Pax6',
    'Pitx2',
    'Prdm1',
    'Pou3f1',
    'Pou5f1',
    'Ripply2',
    'Shh',
    'Six1',
    'Six3',
    'Snai1',
    'Snai2',
    'Sox1',
    'Sox2',
    'Sox3',
    'Sox9',
    'Sox17',
    'T',
    'Tfap2a',
    'Tfap2c',
    'Tbx1',
    'Tbx4',
    'Tbx5',
    'Tbx6',
    'Tbx18',
    'Tcf15',
    'Tcf21',
    'Twist1',
    'Twist2',
    'Uncx',
    'Wnt1',
    'Wnt3',
    'Wnt3a'];
   
    // Ordena el Array Alfabeticamente, es muy facil ;)):
    //array.sort();
   
    addOptions("colorCode", array);
   }
   
// Rutina para agregar opciones a un <select>
function addOptions(domElement, array) {
    var select = document.getElementsByName(domElement)[0];
   
    for (value in array) {
     var option = document.createElement("option");
     option.text = array[value];
     select.add(option);
    }
   }