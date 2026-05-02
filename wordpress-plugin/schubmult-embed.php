<?php
/**
 * Plugin Name: Schubmult Embed
 * Description: Provides a [schubmult] shortcode that embeds the schubmult web
 *              widget (hosted separately) via an iframe. Configure the URL
 *              under Settings → Schubmult Embed.
 * Version:     1.0.0
 * Author:      schubmult
 * License:     MIT
 */

if (!defined('ABSPATH')) { exit; }

const SCHUBMULT_OPTION = 'schubmult_embed_url';
const SCHUBMULT_DEFAULT_URL = 'https://example.pythonanywhere.com/embed';

/**
 * Shortcode: [schubmult height="640" width="100%"]
 *
 * Renders an <iframe> pointing at the configured embed URL.
 */
function schubmult_embed_shortcode($atts) {
    $atts = shortcode_atts(array(
        'height' => 'auto',  // 'auto' = grow to fit content via postMessage
        'width'  => '100%',
        'flavor' => '',     // optional: py | double | q | q_double; appended as ?flavor=...
    ), $atts, 'schubmult');

    $url = trim((string) get_option(SCHUBMULT_OPTION, SCHUBMULT_DEFAULT_URL));
    if (empty($url)) {
        return '<em>Schubmult embed URL is not configured. See Settings → Schubmult Embed.</em>';
    }

    $src = esc_url($url);
    if (!empty($atts['flavor']) && in_array($atts['flavor'], array('py', 'double', 'q', 'q_double'), true)) {
        $sep = (strpos($src, '?') === false) ? '?' : '&';
        $src .= $sep . 'flavor=' . rawurlencode($atts['flavor']);
    }

    $height_attr = strtolower(trim((string) $atts['height']));
    $auto = ($height_attr === 'auto' || $height_attr === '');
    if ($auto) {
        // Initial height before the iframe reports its real size. Kept small
        // so there's no large empty gap on first paint.
        $height = '120';
    } else {
        $height = preg_replace('/[^0-9]/', '', $height_attr);
        if ($height === '') { $height = '640'; }
    }
    $width  = esc_attr($atts['width']);

    // Compute the embed origin so we can validate postMessage senders.
    $origin = '';
    $parts  = wp_parse_url($url);
    if (!empty($parts['scheme']) && !empty($parts['host'])) {
        $origin = $parts['scheme'] . '://' . $parts['host'];
        if (!empty($parts['port'])) { $origin .= ':' . $parts['port']; }
    }

    // Unique id so multiple embeds on one page each resize independently.
    static $counter = 0;
    $counter++;
    $id = 'schubmult-embed-' . $counter;

    $iframe = sprintf(
        '<iframe id="%s" src="%s" width="%s" height="%s" loading="lazy" '
        . 'style="border: 1px solid #ddd; border-radius: 6px; display: block;" '
        . 'sandbox="allow-scripts allow-same-origin allow-forms"></iframe>',
        esc_attr($id), $src, $width, esc_attr($height)
    );

    if (!$auto) {
        return $iframe;
    }

    // Inline auto-resize listener. Scoped to this iframe id and origin.
    $script = sprintf(
        '<script>(function(){'
        . 'var f=document.getElementById(%s);'
        . 'var origin=%s;'
        . 'if(!f)return;'
        . 'window.addEventListener("message",function(ev){'
        .   'if(origin&&ev.origin!==origin)return;'
        .   'var d=ev.data;'
        .   'if(!d||d.type!=="schubmult-resize")return;'
        .   'var h=parseInt(d.height,10);'
        .   'if(h>0&&h<5000)f.style.height=h+"px";'
        . '});'
        . '})();</script>',
        wp_json_encode($id),
        wp_json_encode($origin)
    );

    return $iframe . $script;
}
add_shortcode('schubmult', 'schubmult_embed_shortcode');

/**
 * Settings page.
 */
function schubmult_embed_admin_menu() {
    add_options_page(
        'Schubmult Embed',
        'Schubmult Embed',
        'manage_options',
        'schubmult-embed',
        'schubmult_embed_settings_page'
    );
}
add_action('admin_menu', 'schubmult_embed_admin_menu');

function schubmult_embed_settings_init() {
    register_setting('schubmult_embed', SCHUBMULT_OPTION, array(
        'type'              => 'string',
        'sanitize_callback' => 'esc_url_raw',
        'default'           => SCHUBMULT_DEFAULT_URL,
    ));
}
add_action('admin_init', 'schubmult_embed_settings_init');

function schubmult_embed_settings_page() {
    if (!current_user_can('manage_options')) { return; }
    $url = (string) get_option(SCHUBMULT_OPTION, SCHUBMULT_DEFAULT_URL);
    ?>
    <div class="wrap">
        <h1>Schubmult Embed</h1>
        <p>Set the URL of your hosted schubmult Flask app's <code>/embed</code> endpoint.
           Then drop <code>[schubmult]</code> into any post or page.</p>
        <form action="options.php" method="post">
            <?php settings_fields('schubmult_embed'); ?>
            <table class="form-table" role="presentation">
                <tr>
                    <th scope="row"><label for="schubmult_embed_url">Embed URL</label></th>
                    <td>
                        <input id="schubmult_embed_url" name="<?php echo esc_attr(SCHUBMULT_OPTION); ?>"
                               type="url" value="<?php echo esc_attr($url); ?>"
                               class="regular-text" required>
                        <p class="description">e.g. <code>https://yourname.pythonanywhere.com/embed</code></p>
                    </td>
                </tr>
            </table>
            <?php submit_button(); ?>
        </form>
        <h2>Shortcode reference</h2>
        <pre>[schubmult]                       (auto-resizes to fit content)
[schubmult height="500"]          (fixed pixel height, no auto-resize)
[schubmult flavor="double"]
[schubmult flavor="q_double"]</pre>
        <p>By default the iframe reports its content height back to the page
           and grows/shrinks to fit. Pass an explicit numeric <code>height</code>
           to disable auto-resize.</p>
        <p>The host serving the embed must include your WordPress site in its
           <code>SCHUBMULT_ALLOWED_ORIGINS</code> environment variable so the
           iframe is allowed.</p>
    </div>
    <?php
}
