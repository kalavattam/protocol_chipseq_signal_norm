
# Accepted section-boundary matrix
The first section contains content.

<br />

## Ordinary section without an anchor
This section contains content.

<br />

<a id="ordinary-with-anchor"></a>
## Ordinary section with an anchor
This section contains content.

<br />

## Contentless parent without anchors
### Direct child without anchors
This child contains content.

<br />

<a id="parent-with-anchor"></a>
## Contentless parent with an anchor
<a id="child-with-anchor"></a>
### Direct child with an anchor
This child contains content.

<br />

## Parent containing content
Real parent content preserves the ordinary separator before its child.

<br />

<a id="child-after-content"></a>
### Child after content
This child contains content.

<br />

## Empty sibling

<br />

<a id="empty-sibling-with-anchor"></a>
## Empty sibling with an anchor

<br />

#### Skipped deeper section after an ordinary boundary
This section contains content.

<br />

## Shallower section
<details>
<summary>Details boundary</summary>
Content.
</details>
<br />

<a id="after-details"></a>
## Section after a details close
The direct break after the close is the complete boundary.
