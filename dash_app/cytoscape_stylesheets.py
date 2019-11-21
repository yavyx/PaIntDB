combined = [
    {
        'selector': 'node',
        'style': {
            'content': 'data(label)',
            'border-width': 1,
            'border-color': 'black',
            'min-zoomed-font-size': 3
        }
    },
    {
        'selector': '[significanceSource = "RNASeq"]',
        'style': {
            'background-color': 'green'
        }
    },
    {
        'selector': '[significanceSource = "TnSeq"]',
        'style': {
            'background-color': 'blue'
        }
    },
    {
        'selector': '[significanceSource = "both"]',
        'style': {
            'background-color': 'red'
        }
    }
]

fold_change = [
    {
        'selector': 'node',
        'style': {
            'content': 'data(label)',
            'border-width': 1,
            'border-color': 'black',
            'background-color': 'mapData(log2FoldChange, -4, 4, red, blue)',
            'min-zoomed-font-size': 3
        }
    }
]
