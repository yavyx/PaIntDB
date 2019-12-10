default = [
    {
        'selector': 'node',
        'style': {
            'border-width': 1,
            'border-color': 'black'
        }
    },
    {
        'selector': ':selected',
        'style': {
            'border-width': 8,
            'content': 'data(label)',
            'border-color': 'black',
            'border-opacity': 1,
            'opacity': 1,
            'font-size': 18,
            'z-index': 9999
        }
    },
    {
        'selector': 'edge',
        'style': {
            'width': '2'
        }
    }
]

fold_change = default + \
    [
        {
            'selector': '[significanceSource = "TnSeq"]',
            'style': {
                'background-color': 'gray'
            }
        },
        {
            'selector': '[log2FoldChange < 0]',
            'style': {
                'background-color': '#459cff'  # blue
            }
        },
        {
            'selector': '[log2FoldChange > 0]',
            'style': {
                'background-color': 'red'
            }
        },
    ]


combined = default + \
    [
        {
            'selector': '[significanceSource = "RNASeq"]',
            'style': {
                'background-color': '#26e81c'  # green
            }
        },
        {
            'selector': '[significanceSource = "TnSeq"]',
            'style': {
                'background-color': '#459cff'  # blue
            }
        },
        {
            'selector': '[significanceSource = "both"]',
            'style': {
                'background-color': 'red'
            }
        }
    ]

