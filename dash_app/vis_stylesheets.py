default = [
    {
        'selector': 'node',
        'style': {
            'border-width': 1,
            'border-color': 'black',
            'width': 'mapData(degree, 1, 15, 10, 30)',
            'height': 'mapData(degree, 1, 15, 10, 30)',
            'padding': '10%',
            'font-size': 22,
        }
    },
    {
        'selector': ':selected',
        'style': {
            'border-width': 4,
            'label': 'data(label)',
            'border-color': 'black',
            'border-opacity': 1,
            'opacity': 1,
            'font-weight': 'bold',
            # 'font-size': 26,
            'z-index': 9999  # bring nodes to front
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
                'background-color': '#a6a6a6'  # gray
            }
        },
        {
            'selector': '[log2FoldChange < 0]',
            'style': {
                'background-color': '#0037ff'  # blue
            }
        },
        {
            'selector': '[log2FoldChange > 0]',
            'style': {
                'background-color': '#ff0000'  # red
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
                'background-color': '#ff931f'  # orange
            }
        },
        {
            'selector': '[significanceSource = "both"]',
            'style': {
                'background-color': '#a01cff'  # purple
            }
        }
    ]


def add_labels(stylesheet):
    """Add node labels to the stylesheet."""
    new_stylesheet = stylesheet + \
        [
            {
                'selector': 'node',
                'style': {
                    'label': 'data(label)'
                }
            }
        ]
    return new_stylesheet



