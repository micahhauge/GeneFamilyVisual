/* Author: Micah Hauge
 * Purpose: crete sample json to be visualized
*/

// create base object to be JSONified
var gfv_json = { exon_list: null, family_list: null};

// create exon list and add it to the base object
// FORMAT: {COLUMN, TEXT, LENGTH, FULLY_UTR, [manual_x_adjustment], [manual_y_adjustment]}
gfv_json.exon_list = [
  { column: 1, text: "69", length: 1,  fullyUTR: 1 },
  { column: 1, text: "588", length: 1, fullyUTR: 1 },
  { column: 1, text: "330", length: 1, fullyUTR: 1 },
  { column: 1, text: "533", length: 1, fullyUTR: 1 },
  { column: 2, text: "809", length: 1, fullyUTR: 0 },
  { column: 2, text: "837", length: 1, fullyUTR: 0 },
  { column: 2, text: "859", length: 1, fullyUTR: 0 },
  { column: 2, text: "986", length: 1, fullyUTR: 0 },
  { column: 2, text: "841", length: 1, fullyUTR: 0 },
  { column: 2, text: "541", length: 1, fullyUTR: 0 },
  { column: 2, text: "697", length: 1, fullyUTR: 0 },
  { column: 0, text: "274", length: 2, fullyUTR: 0, manual_x_adjustment: 260, manual_y_adjustment: p.center - p.yPadding * 2 },
  { column: 3, text: "143", length: 1, fullyUTR: 0 },
  { column: 3, text: "142,139", length: 1, fullyUTR: 0 },
  { column: 4, text: "140", length: 1, fullyUTR: 0 },
  { column: 4, text: "132,135", length: 1, fullyUTR: 0 },
  { column: 5, text: "121,127", length: 1, fullyUTR: 0 },
  { column: 6, text: "155", length: 1, fullyUTR: 0 },
  { column: 7, text: "114", length: 1, fullyUTR: 0 },
  { column: 8, text: "132,135,138,141,144", length: 1, fullyUTR: 0 },
  { column: 9, text: "66", length: 1, fullyUTR: 0 },
  { column: 0, text: "132", length: 2, fullyUTR: 0, manual_x_adjustment: 810, manual_y_adjustment: p.center + p.yPadding * .3 },
  { column: 0, text: "72,69", length: 1, fullyUTR: 0, manual_x_adjustment: 950, manual_y_adjustment: p.center + p.yPadding * .5 },
  { column: 10, text: "66", length: 1, fullyUTR: 0 },
  { column: 11, text: "152,155,161,173,176,188", length: 1, fullyUTR: 0 },
  { column: 12, text: "77,80,89", length: 1, fullyUTR: 0 },
  { column: 13, text: "134,137,140", length: 1, fullyUTR: 0 },
  { column: 0, text: "184", length: 2, fullyUTR: 0, manual_x_adjustment: 1165, manual_y_adjustment: p.center + p.yPadding * .3 },
  { column: 14, text: "56", length: 1, fullyUTR: 0 },
  { column: 15, text: "148,151,154,157", length: 1, fullyUTR: 0 },
  { column: 16, text: "74,92,95,98", length: 1, fullyUTR: 0 },
  { column: 0, text: "185", length: 2, fullyUTR: 0, manual_x_adjustment: 1450, manual_y_adjustment: p.center + p.yPadding * .3 },
  { column: 17, text: "111,105,108", length: 1, fullyUTR: 0 },
  { column: 18, text: "162,168,171", length: 1, fullyUTR: 0 },
  { column: 19, text: "141,138", length: 1, fullyUTR: 0 },
  { column: 20, text: "76,79,82", length: 1, fullyUTR: 0 },
  { column: 21, text: "150?", length: 1, fullyUTR: 0 },
  { column: 21, text: "607", length: 1, fullyUTR: 0 },
  { column: 21, text: "150?", length: 1, fullyUTR: 0 },
  { column: 21, text: "150?", length: 1, fullyUTR: 0 },
  { column: 21, text: "150?", length: 1, fullyUTR: 0 },
  { column: 21, text: "150?", length: 1, fullyUTR: 0 },
  { column: 21, text: "2233", length: 1, fullyUTR: 0 },
];

// create a set of colors for family identification
colors = create_random_colors(8);

// create family list
// FORMAT: { EXON_ARRAY, COLOR, FAMILY_NAME }
// gfv_json.family_list = [
//   { exons_in_family: [4, 11, 16, 17, 18, 19, 20, 23, 24, 25, 26, 28, 29, 30, 32, 33, 34, 35, 36], color: colors[0], family_name: "TRHDE"},
//   { exons_in_family: [0, 5, 12, 14, 16, 17, 18, 19, 20, 23, 24, 25, 26, 28, 29, 30, 32, 33, 34, 35, 37], color: colors[1], family_name: "ANPEP"},
//   { exons_in_family: [6, 12, 14, 16, 17, 18, 19, 20, 23, 24, 25, 26, 28, 29, 30, 32, 33, 34, 35, 38], color: colors[2], family_name: "APQ"},
//   { exons_in_family: [7, 13, 15, 16, 17, 18, 19, 20, 23, 24, 25, 26, 28, 29, 30, 32, 33, 34, 35, 39], color: colors[3], family_name: "ENPEP"},
//   { exons_in_family: [3, 8, 13, 15, 16, 17, 18, 19, 21, 24, 25, 27, 29, 31, 33, 34, 35, 40], color: colors[4], family_name: "LNPEP"},
//   { exons_in_family: [2, 9, 13, 15, 16, 17, 18, 19, 21, 22, 24, 25, 27, 29, 31, 33, 34, 35, 41], color: colors[5], family_name: "ERAP1"},
//   { exons_in_family: [1, 10, 13, 15, 16, 17, 18, 19, 21, 22, 24, 25, 27, 29, 31, 33, 34, 35, 42], color: colors[6], family_name: "ERAP2"},
// ]

// NO COLORS:
// create family list
// FORMAT: { EXON_ARRAY, COLOR, FAMILY_NAME }
gfv_json.family_list = [
  { exons_in_family: [4, 11, 16, 17, 18, 19, 20, 23, 24, 25, 26, 28, 29, 30, 32, 33, 34, 35, 36], family_name: "TRHDE"},
  { exons_in_family: [0, 5, 12, 14, 16, 17, 18, 19, 20, 23, 24, 25, 26, 28, 29, 30, 32, 33, 34, 35, 37], family_name: "ANPEP"},
  { exons_in_family: [6, 12, 14, 16, 17, 18, 19, 20, 23, 24, 25, 26, 28, 29, 30, 32, 33, 34, 35, 38], family_name: "APQ"},
  { exons_in_family: [7, 13, 15, 16, 17, 18, 19, 20, 23, 24, 25, 26, 28, 29, 30, 32, 33, 34, 35, 39], family_name: "ENPEP"},
  { exons_in_family: [3, 8, 13, 15, 16, 17, 18, 19, 21, 24, 25, 27, 29, 31, 33, 34, 35, 40], family_name: "LNPEP"},
  { exons_in_family: [2, 9, 13, 15, 16, 17, 18, 19, 21, 22, 24, 25, 27, 29, 31, 33, 34, 35, 41], family_name: "ERAP1"},
  { exons_in_family: [1, 10, 13, 15, 16, 17, 18, 19, 21, 22, 24, 25, 27, 29, 31, 33, 34, 35, 42], family_name: "ERAP2"},
]

// stringify the javascript object
var gfv_json_string = JSON.stringify(gfv_json, null);
